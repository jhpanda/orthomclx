#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  uint32_t query_id;
  uint32_t subject_id;
  uint32_t query_taxon_id;
  uint32_t subject_taxon_id;
  float evalue_mant;
  int32_t evalue_exp;
  float percent_identity;
  float percent_match;
} SimilarityRecordBin;

typedef struct {
  uint32_t record_index;
} RecordRef;

typedef struct {
  uint32_t seq_a;
  uint32_t seq_b;
  uint32_t taxon_a;
  uint32_t taxon_b;
  double unnormalized_score;
  double normalized_score;
} CoorthologEdge;

typedef struct {
  SimilarityRecordBin *items;
  size_t len;
  size_t cap;
} RecordVec;

typedef struct {
  RecordRef *items;
  size_t len;
  size_t cap;
} RefVec;

typedef struct {
  CoorthologEdge *items;
  size_t len;
  size_t cap;
} EdgeVec;

typedef struct {
  uint32_t taxon_a;
  uint32_t taxon_b;
  double sum;
  uint32_t count;
} TaxonPairAgg;

typedef struct {
  TaxonPairAgg *items;
  size_t len;
  size_t cap;
} TaxonPairAggVec;

typedef struct {
  char **items;
  size_t len;
  size_t cap;
} StringVec;

typedef struct {
  char *name;
  uint32_t index;
} NameIndex;

typedef struct {
  NameIndex *items;
  size_t len;
  size_t cap;
} NameIndexVec;

typedef struct {
  uint32_t *items;
  size_t len;
  size_t cap;
} UIntVec;

typedef struct {
  UIntVec *items;
  size_t len;
} NeighborGraph;

typedef struct {
  uint64_t *items;
  size_t len;
  size_t cap;
} PairSet;

static RecordVec g_records = {0};

static void die(const char *message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
}

static void *xmalloc(size_t size) {
  void *ptr = malloc(size);
  if (!ptr) die("Out of memory");
  return ptr;
}

static void *xrealloc(void *ptr, size_t size) {
  void *out = realloc(ptr, size);
  if (!out) die("Out of memory");
  return out;
}

static char *xstrdup(const char *src) {
  size_t len = strlen(src);
  char *dst = (char *)xmalloc(len + 1);
  memcpy(dst, src, len + 1);
  return dst;
}

static int read_line_alloc(FILE *handle, char **buffer, size_t *capacity) {
  if (*buffer == NULL || *capacity == 0) {
    *capacity = 1024;
    *buffer = (char *)xmalloc(*capacity);
  }
  size_t length = 0;
  for (;;) {
    int ch = fgetc(handle);
    if (ch == EOF) {
      if (length == 0) return 0;
      break;
    }
    if (length + 1 >= *capacity) {
      *capacity *= 2;
      *buffer = (char *)xrealloc(*buffer, *capacity);
    }
    (*buffer)[length++] = (char)ch;
    if (ch == '\n') break;
  }
  (*buffer)[length] = '\0';
  return 1;
}

static void split_first_two_tab_fields(char *line, char **first, char **second) {
  char *tab = strchr(line, '\t');
  if (!tab) {
    *first = NULL;
    *second = NULL;
    return;
  }
  *tab = '\0';
  *first = line;
  *second = tab + 1;
  char *end = *second;
  while (*end && *end != '\t' && *end != '\n' && *end != '\r') {
    end++;
  }
  *end = '\0';
}

static void push_record(RecordVec *vec, SimilarityRecordBin item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 1024;
    vec->items = (SimilarityRecordBin *)xrealloc(vec->items, vec->cap * sizeof(SimilarityRecordBin));
  }
  vec->items[vec->len++] = item;
}

static void push_edge(EdgeVec *vec, CoorthologEdge item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 256;
    vec->items = (CoorthologEdge *)xrealloc(vec->items, vec->cap * sizeof(CoorthologEdge));
  }
  vec->items[vec->len++] = item;
}

static void push_taxon_pair_agg(TaxonPairAggVec *vec, TaxonPairAgg item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 32;
    vec->items = (TaxonPairAgg *)xrealloc(vec->items, vec->cap * sizeof(TaxonPairAgg));
  }
  vec->items[vec->len++] = item;
}

static void push_uint(UIntVec *vec, uint32_t value) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 8;
    vec->items = (uint32_t *)xrealloc(vec->items, vec->cap * sizeof(uint32_t));
  }
  vec->items[vec->len++] = value;
}

static uint64_t pair_key(uint32_t seq_a, uint32_t seq_b) {
  uint32_t first = seq_a < seq_b ? seq_a : seq_b;
  uint32_t second = seq_a < seq_b ? seq_b : seq_a;
  return ((uint64_t)first << 32) | (uint64_t)second;
}

static size_t next_power_of_two(size_t value) {
  size_t size = 1;
  while (size < value) size <<= 1;
  return size;
}

static void pair_set_init(PairSet *set, size_t expected_items) {
  set->cap = next_power_of_two(expected_items ? expected_items * 4 : 1024);
  set->len = 0;
  set->items = (uint64_t *)calloc(set->cap, sizeof(uint64_t));
  if (!set->items) die("Out of memory");
}

static bool pair_set_contains(const PairSet *set, uint64_t key) {
  size_t mask = set->cap - 1;
  size_t slot = (size_t)((key * 11400714819323198485ull) & mask);
  while (set->items[slot] != 0) {
    if (set->items[slot] == key) return true;
    slot = (slot + 1) & mask;
  }
  return false;
}

static void pair_set_grow(PairSet *set) {
  PairSet grown = {0};
  pair_set_init(&grown, set->cap);
  for (size_t i = 0; i < set->cap; i++) {
    uint64_t key = set->items[i];
    if (!key) continue;
    size_t mask = grown.cap - 1;
    size_t slot = (size_t)((key * 11400714819323198485ull) & mask);
    while (grown.items[slot] != 0) slot = (slot + 1) & mask;
    grown.items[slot] = key;
    grown.len += 1;
  }
  free(set->items);
  *set = grown;
}

static void pair_set_add(PairSet *set, uint64_t key) {
  if ((set->len + 1) * 10 >= set->cap * 7) pair_set_grow(set);
  size_t mask = set->cap - 1;
  size_t slot = (size_t)((key * 11400714819323198485ull) & mask);
  while (set->items[slot] != 0) {
    if (set->items[slot] == key) return;
    slot = (slot + 1) & mask;
  }
  set->items[slot] = key;
  set->len += 1;
}

static int cmp_name_index(const void *left, const void *right) {
  const NameIndex *a = (const NameIndex *)left;
  const NameIndex *b = (const NameIndex *)right;
  return strcmp(a->name, b->name);
}

static int cmp_ref_query_subject(const void *left, const void *right) {
  const SimilarityRecordBin *a = &g_records.items[((const RecordRef *)left)->record_index];
  const SimilarityRecordBin *b = &g_records.items[((const RecordRef *)right)->record_index];
  if (a->query_id != b->query_id) return a->query_id < b->query_id ? -1 : 1;
  if (a->subject_id != b->subject_id) return a->subject_id < b->subject_id ? -1 : 1;
  return 0;
}

static int cmp_edge_output(const void *left, const void *right) {
  const CoorthologEdge *a = (const CoorthologEdge *)left;
  const CoorthologEdge *b = (const CoorthologEdge *)right;
  if (a->taxon_a != b->taxon_a) return a->taxon_a < b->taxon_a ? -1 : 1;
  if (a->taxon_b != b->taxon_b) return a->taxon_b < b->taxon_b ? -1 : 1;
  if (a->seq_a != b->seq_a) return a->seq_a < b->seq_a ? -1 : 1;
  if (a->seq_b != b->seq_b) return a->seq_b < b->seq_b ? -1 : 1;
  return 0;
}

static RecordVec read_binary_records(const char *path) {
  RecordVec records = {0};
  FILE *handle = fopen(path, "rb");
  if (!handle) die("Could not open similarities.bin");
  SimilarityRecordBin record;
  int32_t min_nonzero_exp = INT32_MAX;
  while (fread(&record, sizeof(SimilarityRecordBin), 1, handle) == 1) {
    if (record.evalue_mant != 0.0f && record.evalue_exp < min_nonzero_exp) {
      min_nonzero_exp = record.evalue_exp;
    }
    push_record(&records, record);
  }
  fclose(handle);
  if (min_nonzero_exp != INT32_MAX) {
    int32_t adjusted_zero_exp = min_nonzero_exp - 1;
    for (size_t i = 0; i < records.len; i++) {
      if (records.items[i].evalue_mant == 0.0f && records.items[i].evalue_exp == 0) {
        records.items[i].evalue_exp = adjusted_zero_exp;
      }
    }
  }
  return records;
}

static StringVec read_index_values(const char *path) {
  StringVec values = {0};
  FILE *handle = fopen(path, "r");
  if (!handle) die("Could not open index file");
  char *line = NULL;
  size_t cap = 0;
  while (read_line_alloc(handle, &line, &cap)) {
    char *tab = strchr(line, '\t');
    if (!tab) die("Invalid index file");
    char *value = tab + 1;
    char *newline = strchr(value, '\n');
    if (newline) *newline = '\0';
    if (values.len == values.cap) {
      values.cap = values.cap ? values.cap * 2 : 256;
      values.items = (char **)xrealloc(values.items, values.cap * sizeof(char *));
    }
    values.items[values.len++] = xstrdup(value);
  }
  free(line);
  fclose(handle);
  return values;
}

static NameIndexVec build_name_index(StringVec *proteins) {
  NameIndexVec map = {0};
  map.len = proteins->len;
  map.cap = proteins->len;
  map.items = (NameIndex *)xmalloc(map.len * sizeof(NameIndex));
  for (size_t i = 0; i < proteins->len; i++) {
    map.items[i].name = proteins->items[i];
    map.items[i].index = (uint32_t)i;
  }
  qsort(map.items, map.len, sizeof(NameIndex), cmp_name_index);
  return map;
}

static int name_index_lookup(NameIndexVec *map, const char *name) {
  size_t left = 0;
  size_t right = map->len;
  while (left < right) {
    size_t mid = left + (right - left) / 2;
    int cmp = strcmp(map->items[mid].name, name);
    if (cmp < 0) left = mid + 1;
    else right = mid;
  }
  if (left < map->len && strcmp(map->items[left].name, name) == 0) return (int)map->items[left].index;
  return -1;
}

static RefVec build_sorted_refs(size_t record_count) {
  RefVec refs = {0};
  refs.items = (RecordRef *)xmalloc(record_count * sizeof(RecordRef));
  refs.len = record_count;
  refs.cap = record_count;
  for (size_t i = 0; i < record_count; i++) refs.items[i].record_index = (uint32_t)i;
  qsort(refs.items, refs.len, sizeof(RecordRef), cmp_ref_query_subject);
  return refs;
}

static SimilarityRecordBin *find_record(RefVec *refs, uint32_t query_id, uint32_t subject_id) {
  size_t left = 0;
  size_t right = refs->len;
  while (left < right) {
    size_t mid = left + (right - left) / 2;
    SimilarityRecordBin *record = &g_records.items[refs->items[mid].record_index];
    if (record->query_id < query_id || (record->query_id == query_id && record->subject_id < subject_id)) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  if (left < refs->len) {
    SimilarityRecordBin *record = &g_records.items[refs->items[left].record_index];
    if (record->query_id == query_id && record->subject_id == subject_id) return record;
  }
  return NULL;
}

static bool passes_evalue_cutoff(const SimilarityRecordBin *record, float cutoff_mant, int32_t cutoff_exp) {
  if (record->evalue_mant == 0.0f) return true;
  if (record->evalue_exp != cutoff_exp) return record->evalue_exp < cutoff_exp;
  return record->evalue_mant <= cutoff_mant;
}

static bool passes_thresholds(const SimilarityRecordBin *record, float percent_match_cutoff, float cutoff_mant, int32_t cutoff_exp) {
  return record->percent_match >= percent_match_cutoff && passes_evalue_cutoff(record, cutoff_mant, cutoff_exp);
}

static double score_from_records(const SimilarityRecordBin *left, const SimilarityRecordBin *right, double zero_cutoff) {
  if (left->evalue_mant < zero_cutoff || right->evalue_mant < zero_cutoff) {
    return ((double)left->evalue_exp + (double)right->evalue_exp) / -2.0;
  }
  return (log10((double)left->evalue_mant * (double)right->evalue_mant) +
          (double)left->evalue_exp + (double)right->evalue_exp) /
         -2.0;
}

static NeighborGraph init_neighbor_graph(size_t node_count) {
  NeighborGraph graph = {0};
  graph.items = (UIntVec *)calloc(node_count, sizeof(UIntVec));
  if (!graph.items) die("Out of memory");
  graph.len = node_count;
  return graph;
}

static void add_neighbor(NeighborGraph *graph, uint32_t node, uint32_t neighbor) {
  UIntVec *entry = &graph->items[node];
  for (size_t i = 0; i < entry->len; i++) {
    if (entry->items[i] == neighbor) return;
  }
  push_uint(entry, neighbor);
}

static UIntVec *find_neighbors(NeighborGraph *graph, uint32_t node) {
  if (node >= graph->len) return NULL;
  return &graph->items[node];
}

static void load_pair_graph(const char *path, NameIndexVec *protein_map, NeighborGraph *graph) {
  FILE *handle = fopen(path, "r");
  if (!handle) die("Could not open edge file");
  char *line = NULL;
  size_t cap = 0;
  while (read_line_alloc(handle, &line, &cap)) {
    char *seq_a = NULL;
    char *seq_b = NULL;
    split_first_two_tab_fields(line, &seq_a, &seq_b);
    if (!seq_a || !seq_b) die("Invalid edge file");
    int a = name_index_lookup(protein_map, seq_a);
    int b = name_index_lookup(protein_map, seq_b);
    if (a < 0 || b < 0) die("Protein from edge file not found in proteins.tsv");
    add_neighbor(graph, (uint32_t)a, (uint32_t)b);
    add_neighbor(graph, (uint32_t)b, (uint32_t)a);
  }
  free(line);
  fclose(handle);
}

static bool edge_in_graph(NeighborGraph *graph, uint32_t seq_a, uint32_t seq_b) {
  UIntVec *neighbors = find_neighbors(graph, seq_a);
  if (!neighbors) return false;
  for (size_t i = 0; i < neighbors->len; i++) {
    if (neighbors->items[i] == seq_b) return true;
  }
  return false;
}

static void maybe_add_coortholog(
    EdgeVec *edges,
    PairSet *seen_pairs,
    NeighborGraph *orth_graph,
    RefVec *by_query_subject,
    StringVec *proteins,
    uint32_t seq_a,
    uint32_t seq_b,
    float percent_match_cutoff,
    float cutoff_mant,
    int32_t cutoff_exp
) {
  uint32_t first = seq_a;
  uint32_t second = seq_b;
  if (strcmp(proteins->items[first], proteins->items[second]) > 0) {
    first = seq_b;
    second = seq_a;
  }
  if (first == second) return;
  if (edge_in_graph(orth_graph, first, second)) return;
  uint64_t key = pair_key(first, second);
  if (pair_set_contains(seen_pairs, key)) return;

  SimilarityRecordBin *left = find_record(by_query_subject, first, second);
  SimilarityRecordBin *right = find_record(by_query_subject, second, first);
  if (!left || !right) return;
  if (!passes_thresholds(left, percent_match_cutoff, cutoff_mant, cutoff_exp)) return;
  if (!passes_thresholds(right, percent_match_cutoff, cutoff_mant, cutoff_exp)) return;

  CoorthologEdge edge;
  edge.seq_a = first;
  edge.seq_b = second;
  edge.taxon_a = left->query_taxon_id;
  edge.taxon_b = left->subject_taxon_id;
  edge.unnormalized_score = score_from_records(left, right, 0.00001);
  edge.normalized_score = 0.0;
  pair_set_add(seen_pairs, key);
  push_edge(edges, edge);
}

static TaxonPairAgg *get_taxon_pair_agg(TaxonPairAggVec *aggs, uint32_t taxon_a, uint32_t taxon_b) {
  uint32_t first = taxon_a < taxon_b ? taxon_a : taxon_b;
  uint32_t second = taxon_a < taxon_b ? taxon_b : taxon_a;
  for (size_t i = 0; i < aggs->len; i++) {
    if (aggs->items[i].taxon_a == first && aggs->items[i].taxon_b == second) return &aggs->items[i];
  }
  TaxonPairAgg agg;
  agg.taxon_a = first;
  agg.taxon_b = second;
  agg.sum = 0.0;
  agg.count = 0;
  push_taxon_pair_agg(aggs, agg);
  return &aggs->items[aggs->len - 1];
}

static EdgeVec build_coortholog_edges(
    NeighborGraph *orth_graph,
    NeighborGraph *in_graph,
    RefVec *by_query_subject,
    StringVec *proteins,
    float percent_match_cutoff,
    float cutoff_mant,
    int32_t cutoff_exp,
    uint32_t shard_index,
    uint32_t shard_count
) {
  EdgeVec edges = {0};
  PairSet seen_pairs = {0};
  pair_set_init(&seen_pairs, 262144);
  uint32_t *target_marks = (uint32_t *)calloc(in_graph->len, sizeof(uint32_t));
  if (!target_marks) die("Out of memory");
  uint32_t mark_generation = 1;
  UIntVec target_nodes = {0};
  uint32_t processed_sources = 0;

  for (uint32_t source = 0; source < in_graph->len; source++) {
    if (shard_count > 1 && (source % shard_count) != shard_index) continue;
    processed_sources += 1;
    if (processed_sources % 10000 == 0) {
      fprintf(stderr, "[indexed_coorthologs] shard %u/%u processed %u source proteins, %zu coorthologs so far\n",
              shard_index + 1, shard_count, processed_sources, edges.len);
    }
    UIntVec *in_neighbors = &in_graph->items[source];
    if (in_neighbors->len == 0) continue;
    UIntVec *orth_neighbors = find_neighbors(orth_graph, source);
    if (!orth_neighbors || orth_neighbors->len == 0) continue;

    if (mark_generation == UINT32_MAX) {
      memset(target_marks, 0, in_graph->len * sizeof(uint32_t));
      mark_generation = 1;
    }
    target_nodes.len = 0;

    for (size_t b = 0; b < orth_neighbors->len; b++) {
      uint32_t orth_neighbor = orth_neighbors->items[b];
      if (target_marks[orth_neighbor] != mark_generation) {
        target_marks[orth_neighbor] = mark_generation;
        push_uint(&target_nodes, orth_neighbor);
      }
      UIntVec *right_in = find_neighbors(in_graph, orth_neighbor);
      if (!right_in || right_in->len == 0) continue;
      for (size_t c = 0; c < right_in->len; c++) {
        uint32_t candidate = right_in->items[c];
        if (target_marks[candidate] != mark_generation) {
          target_marks[candidate] = mark_generation;
          push_uint(&target_nodes, candidate);
        }
      }
    }
    mark_generation += 1;

    for (size_t a = 0; a < in_neighbors->len; a++) {
      for (size_t t = 0; t < target_nodes.len; t++) {
        maybe_add_coortholog(
            &edges, &seen_pairs, orth_graph, by_query_subject, proteins,
            in_neighbors->items[a], target_nodes.items[t],
            percent_match_cutoff, cutoff_mant, cutoff_exp);
      }
    }
  }

  TaxonPairAggVec aggs = {0};
  for (size_t i = 0; i < edges.len; i++) {
    TaxonPairAgg *agg = get_taxon_pair_agg(&aggs, edges.items[i].taxon_a, edges.items[i].taxon_b);
    agg->sum += edges.items[i].unnormalized_score;
    agg->count += 1;
  }
  for (size_t i = 0; i < edges.len; i++) {
    TaxonPairAgg *agg = get_taxon_pair_agg(&aggs, edges.items[i].taxon_a, edges.items[i].taxon_b);
    edges.items[i].normalized_score = edges.items[i].unnormalized_score / (agg->sum / (double)agg->count);
  }
  free(aggs.items);
  free(seen_pairs.items);
  free(target_marks);
  free(target_nodes.items);

  qsort(edges.items, edges.len, sizeof(CoorthologEdge), cmp_edge_output);
  return edges;
}

static double round_score(double score) {
  return floor(score * 1000.0 + 0.5) / 1000.0;
}

static void make_dir(const char *path) {
  char command[4096];
  snprintf(command, sizeof(command), "mkdir -p \"%s\"", path);
  if (system(command) != 0) die("Could not create output directory");
}

static void write_coorthologs(const char *path, EdgeVec *edges, StringVec *proteins) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open coortholog output");
  for (size_t i = 0; i < edges->len; i++) {
    fprintf(handle, "%s\t%s\t%.3f\n",
            proteins->items[edges->items[i].seq_a],
            proteins->items[edges->items[i].seq_b],
            round_score(edges->items[i].normalized_score));
  }
  fclose(handle);
}

static void write_coorthologs_raw(const char *path, EdgeVec *edges, StringVec *proteins, StringVec *taxa) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open coortholog raw output");
  for (size_t i = 0; i < edges->len; i++) {
    fprintf(handle, "%s\t%s\t%s\t%s\t%.12f\n",
            proteins->items[edges->items[i].seq_a],
            proteins->items[edges->items[i].seq_b],
            taxa->items[edges->items[i].taxon_a],
            taxa->items[edges->items[i].taxon_b],
            edges->items[i].unnormalized_score);
  }
  fclose(handle);
}

static void write_summary(FILE *handle, size_t record_count, size_t orth_nodes, size_t in_nodes, size_t edge_count, size_t protein_count, size_t taxon_count) {
  fprintf(handle, "records\t%zu\n", record_count);
  fprintf(handle, "ortholog_nodes\t%zu\n", orth_nodes);
  fprintf(handle, "inparalog_nodes\t%zu\n", in_nodes);
  fprintf(handle, "coortholog_edges\t%zu\n", edge_count);
  fprintf(handle, "proteins\t%zu\n", protein_count);
  fprintf(handle, "taxa\t%zu\n", taxon_count);
}

static void free_strings(StringVec *values) {
  for (size_t i = 0; i < values->len; i++) free(values->items[i]);
  free(values->items);
}

static void free_graph(NeighborGraph *graph) {
  for (size_t i = 0; i < graph->len; i++) {
    free(graph->items[i].items);
  }
  free(graph->items);
}

static void usage(void) {
  fprintf(stderr, "usage: indexed_coorthologs similarities.bin proteins.tsv taxa.tsv orthologs.txt inparalogs.txt out_dir percent_match_cutoff evalue_cutoff_mant evalue_cutoff_exp [shard_index shard_count]\n");
}

int main(int argc, char **argv) {
  if (argc != 10 && argc != 12) {
    usage();
    return 1;
  }

  const char *binary_path = argv[1];
  const char *proteins_path = argv[2];
  const char *taxa_path = argv[3];
  const char *orthologs_path = argv[4];
  const char *inparalogs_path = argv[5];
  const char *out_dir = argv[6];
  float percent_match_cutoff = (float)atof(argv[7]);
  float cutoff_mant = (float)atof(argv[8]);
  int32_t cutoff_exp = (int32_t)atoi(argv[9]);
  uint32_t shard_index = 0;
  uint32_t shard_count = 1;
  if (argc == 12) {
    shard_index = (uint32_t)atoi(argv[10]);
    shard_count = (uint32_t)atoi(argv[11]);
    if (shard_count == 0 || shard_index >= shard_count) die("Invalid shard arguments");
  }

  g_records = read_binary_records(binary_path);
  StringVec proteins = read_index_values(proteins_path);
  StringVec taxa = read_index_values(taxa_path);
  NameIndexVec protein_map = build_name_index(&proteins);
  RefVec by_query_subject = build_sorted_refs(g_records.len);
  NeighborGraph orth_graph = init_neighbor_graph(proteins.len);
  NeighborGraph in_graph = init_neighbor_graph(proteins.len);
  load_pair_graph(orthologs_path, &protein_map, &orth_graph);
  load_pair_graph(inparalogs_path, &protein_map, &in_graph);
  EdgeVec edges = build_coortholog_edges(&orth_graph, &in_graph, &by_query_subject, &proteins, percent_match_cutoff, cutoff_mant, cutoff_exp, shard_index, shard_count);

  make_dir(out_dir);
  char out_path[4096];
  char raw_path[4096];
  char summary_path[4096];
  snprintf(out_path, sizeof(out_path), "%s/coorthologs.indexed.txt", out_dir);
  snprintf(raw_path, sizeof(raw_path), "%s/coorthologs.indexed.raw.tsv", out_dir);
  snprintf(summary_path, sizeof(summary_path), "%s/coorthologs.indexed.summary.tsv", out_dir);
  write_coorthologs(out_path, &edges, &proteins);
  write_coorthologs_raw(raw_path, &edges, &proteins, &taxa);
  FILE *summary = fopen(summary_path, "w");
  if (!summary) die("Could not open summary output");
  write_summary(summary, g_records.len, orth_graph.len, in_graph.len, edges.len, proteins.len, taxa.len);
  fclose(summary);

  free(edges.items);
  free_graph(&orth_graph);
  free_graph(&in_graph);
  free(protein_map.items);
  free(by_query_subject.items);
  free(g_records.items);
  free_strings(&proteins);
  free_strings(&taxa);
  return 0;
}
