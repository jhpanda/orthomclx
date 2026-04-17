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
  uint32_t query_id;
  uint32_t start;
  uint32_t end;
  uint32_t best_inter_record_index;
} QueryGroup;

typedef struct {
  uint32_t seq_a;
  uint32_t seq_b;
  uint32_t taxon_id;
  double unnormalized_score;
  double normalized_score;
} InparalogEdge;

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
  QueryGroup *items;
  size_t len;
  size_t cap;
} GroupVec;

typedef struct {
  InparalogEdge *items;
  size_t len;
  size_t cap;
} EdgeVec;

typedef struct {
  uint32_t taxon_id;
  double all_sum;
  uint32_t all_count;
  double orth_sum;
  uint32_t orth_count;
} TaxonAgg;

typedef struct {
  TaxonAgg *items;
  size_t len;
  size_t cap;
} TaxonAggVec;

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

static void push_record(RecordVec *vec, SimilarityRecordBin item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 1024;
    vec->items = (SimilarityRecordBin *)xrealloc(vec->items, vec->cap * sizeof(SimilarityRecordBin));
  }
  vec->items[vec->len++] = item;
}

static void push_group(GroupVec *vec, QueryGroup item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 256;
    vec->items = (QueryGroup *)xrealloc(vec->items, vec->cap * sizeof(QueryGroup));
  }
  vec->items[vec->len++] = item;
}

static void push_edge(EdgeVec *vec, InparalogEdge item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 256;
    vec->items = (InparalogEdge *)xrealloc(vec->items, vec->cap * sizeof(InparalogEdge));
  }
  vec->items[vec->len++] = item;
}

static void push_taxon_agg(TaxonAggVec *vec, TaxonAgg item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 16;
    vec->items = (TaxonAgg *)xrealloc(vec->items, vec->cap * sizeof(TaxonAgg));
  }
  vec->items[vec->len++] = item;
}

static int cmp_name_index(const void *left, const void *right) {
  const NameIndex *a = (const NameIndex *)left;
  const NameIndex *b = (const NameIndex *)right;
  return strcmp(a->name, b->name);
}

static int cmp_ref_query_evalue(const void *left, const void *right) {
  const SimilarityRecordBin *a = &g_records.items[((const RecordRef *)left)->record_index];
  const SimilarityRecordBin *b = &g_records.items[((const RecordRef *)right)->record_index];
  if (a->query_id != b->query_id) return a->query_id < b->query_id ? -1 : 1;
  if (a->evalue_exp != b->evalue_exp) return a->evalue_exp < b->evalue_exp ? -1 : 1;
  if (a->evalue_mant < b->evalue_mant) return -1;
  if (a->evalue_mant > b->evalue_mant) return 1;
  if (a->subject_id != b->subject_id) return a->subject_id < b->subject_id ? -1 : 1;
  return 0;
}

static int cmp_ref_query_subject(const void *left, const void *right) {
  const SimilarityRecordBin *a = &g_records.items[((const RecordRef *)left)->record_index];
  const SimilarityRecordBin *b = &g_records.items[((const RecordRef *)right)->record_index];
  if (a->query_id != b->query_id) return a->query_id < b->query_id ? -1 : 1;
  if (a->subject_id != b->subject_id) return a->subject_id < b->subject_id ? -1 : 1;
  return 0;
}

static int cmp_edge_output(const void *left, const void *right) {
  const InparalogEdge *a = (const InparalogEdge *)left;
  const InparalogEdge *b = (const InparalogEdge *)right;
  if (a->taxon_id != b->taxon_id) return a->taxon_id < b->taxon_id ? -1 : 1;
  if (a->seq_a != b->seq_a) return a->seq_a < b->seq_a ? -1 : 1;
  if (a->seq_b != b->seq_b) return a->seq_b < b->seq_b ? -1 : 1;
  return 0;
}

static bool passes_thresholds(const SimilarityRecordBin *record, float percent_match_cutoff, int32_t evalue_exp_cutoff) {
  return record->percent_match >= percent_match_cutoff && record->evalue_exp <= evalue_exp_cutoff;
}

static double score_from_records(const SimilarityRecordBin *left, const SimilarityRecordBin *right, double zero_cutoff) {
  if (left->evalue_mant < zero_cutoff || right->evalue_mant < zero_cutoff) {
    return ((double)left->evalue_exp + (double)right->evalue_exp) / -2.0;
  }
  return (log10((double)left->evalue_mant * (double)right->evalue_mant) +
          (double)left->evalue_exp + (double)right->evalue_exp) /
         -2.0;
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
  while (getline(&line, &cap, handle) != -1) {
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

static RefVec build_sorted_refs(size_t record_count, int (*cmp)(const void *, const void *)) {
  RefVec refs = {0};
  refs.items = (RecordRef *)xmalloc(record_count * sizeof(RecordRef));
  refs.len = record_count;
  refs.cap = record_count;
  for (size_t i = 0; i < record_count; i++) refs.items[i].record_index = (uint32_t)i;
  qsort(refs.items, refs.len, sizeof(RecordRef), cmp);
  return refs;
}

static GroupVec build_query_groups(RefVec *by_query_evalue) {
  GroupVec groups = {0};
  size_t i = 0;
  while (i < by_query_evalue->len) {
    const SimilarityRecordBin *first = &g_records.items[by_query_evalue->items[i].record_index];
    uint32_t query_id = first->query_id;
    size_t start = i;
    size_t end = i;
    uint32_t best_inter_record_index = UINT32_MAX;
    while (end < by_query_evalue->len) {
      const SimilarityRecordBin *record = &g_records.items[by_query_evalue->items[end].record_index];
      if (record->query_id != query_id) break;
      if (record->query_taxon_id != record->subject_taxon_id && best_inter_record_index == UINT32_MAX) {
        best_inter_record_index = by_query_evalue->items[end].record_index;
      }
      end++;
    }
    QueryGroup group;
    group.query_id = query_id;
    group.start = (uint32_t)start;
    group.end = (uint32_t)end;
    group.best_inter_record_index = best_inter_record_index;
    push_group(&groups, group);
    i = end;
  }
  return groups;
}

static QueryGroup *find_query_group(GroupVec *groups, uint32_t query_id) {
  size_t left = 0;
  size_t right = groups->len;
  while (left < right) {
    size_t mid = left + (right - left) / 2;
    if (groups->items[mid].query_id < query_id) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  if (left < groups->len && groups->items[left].query_id == query_id) return &groups->items[left];
  return NULL;
}

static SimilarityRecordBin *find_record_by_query_subject(RefVec *refs, uint32_t query_id, uint32_t subject_id) {
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

static bool keep_intra_record(const SimilarityRecordBin *record, QueryGroup *group, float percent_match_cutoff, int32_t evalue_exp_cutoff) {
  if (record->query_taxon_id != record->subject_taxon_id) return false;
  if (record->query_id == record->subject_id) return false;
  if (!passes_thresholds(record, percent_match_cutoff, evalue_exp_cutoff)) return false;
  if (!group || group->best_inter_record_index == UINT32_MAX) return true;

  const SimilarityRecordBin *best_inter = &g_records.items[group->best_inter_record_index];
  if (record->evalue_mant < 0.001f) return true;
  if (record->evalue_exp < best_inter->evalue_exp) return true;
  if (record->evalue_exp == best_inter->evalue_exp && record->evalue_mant <= best_inter->evalue_mant) return true;
  return false;
}

static TaxonAgg *get_taxon_agg(TaxonAggVec *aggs, uint32_t taxon_id) {
  for (size_t i = 0; i < aggs->len; i++) {
    if (aggs->items[i].taxon_id == taxon_id) return &aggs->items[i];
  }
  TaxonAgg agg;
  agg.taxon_id = taxon_id;
  agg.all_sum = 0.0;
  agg.all_count = 0;
  agg.orth_sum = 0.0;
  agg.orth_count = 0;
  push_taxon_agg(aggs, agg);
  return &aggs->items[aggs->len - 1];
}

static bool *load_ortholog_id_mask(const char *orthologs_path, NameIndexVec *protein_map, size_t protein_count) {
  bool *mask = (bool *)calloc(protein_count, sizeof(bool));
  if (!mask) die("Out of memory");
  FILE *handle = fopen(orthologs_path, "r");
  if (!handle) die("Could not open orthologs file");
  char *line = NULL;
  size_t cap = 0;
  while (getline(&line, &cap, handle) != -1) {
    char *saveptr = NULL;
    char *seq_a = strtok_r(line, "\t\n", &saveptr);
    char *seq_b = strtok_r(NULL, "\t\n", &saveptr);
    if (!seq_a || !seq_b) die("Invalid ortholog file");
    int a = name_index_lookup(protein_map, seq_a);
    int b = name_index_lookup(protein_map, seq_b);
    if (a < 0 || b < 0) die("Ortholog protein not found in proteins.tsv");
    if ((size_t)a >= protein_count || (size_t)b >= protein_count) die("Invalid protein index");
    mask[a] = true;
    mask[b] = true;
  }
  free(line);
  fclose(handle);
  return mask;
}

static EdgeVec build_inparalog_edges(GroupVec *groups, RefVec *by_query_subject, StringVec *proteins, bool *ortholog_mask, float percent_match_cutoff, int32_t evalue_exp_cutoff, uint32_t shard_index, uint32_t shard_count) {
  EdgeVec edges = {0};
  for (size_t i = 0; i < g_records.len; i++) {
    SimilarityRecordBin *left = &g_records.items[i];
    if (shard_count > 1 && (left->query_id % shard_count) != shard_index) continue;
    QueryGroup *left_group = find_query_group(groups, left->query_id);
    if (!keep_intra_record(left, left_group, percent_match_cutoff, evalue_exp_cutoff)) continue;
    if (strcmp(proteins->items[left->query_id], proteins->items[left->subject_id]) >= 0) continue;

    SimilarityRecordBin *right = find_record_by_query_subject(by_query_subject, left->subject_id, left->query_id);
    if (!right) continue;
    QueryGroup *right_group = find_query_group(groups, right->query_id);
    if (!keep_intra_record(right, right_group, percent_match_cutoff, evalue_exp_cutoff)) continue;

    InparalogEdge edge;
    if (strcmp(proteins->items[left->query_id], proteins->items[left->subject_id]) < 0) {
      edge.seq_a = left->query_id;
      edge.seq_b = left->subject_id;
    } else {
      edge.seq_a = left->subject_id;
      edge.seq_b = left->query_id;
    }
    edge.taxon_id = left->query_taxon_id;
    edge.unnormalized_score = score_from_records(left, right, 0.01);
    edge.normalized_score = 0.0;
    push_edge(&edges, edge);
  }

  TaxonAggVec aggs = {0};
  for (size_t i = 0; i < edges.len; i++) {
    TaxonAgg *agg = get_taxon_agg(&aggs, edges.items[i].taxon_id);
    agg->all_sum += edges.items[i].unnormalized_score;
    agg->all_count += 1;
    if (ortholog_mask[edges.items[i].seq_a] || ortholog_mask[edges.items[i].seq_b]) {
      agg->orth_sum += edges.items[i].unnormalized_score;
      agg->orth_count += 1;
    }
  }
  for (size_t i = 0; i < edges.len; i++) {
    TaxonAgg *agg = get_taxon_agg(&aggs, edges.items[i].taxon_id);
    double avg = agg->orth_count ? (agg->orth_sum / (double)agg->orth_count) : (agg->all_sum / (double)agg->all_count);
    edges.items[i].normalized_score = edges.items[i].unnormalized_score / avg;
  }
  free(aggs.items);

  qsort(edges.items, edges.len, sizeof(InparalogEdge), cmp_edge_output);
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

static void write_inparalogs(const char *path, EdgeVec *edges, StringVec *proteins) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open inparalog output");
  for (size_t i = 0; i < edges->len; i++) {
    fprintf(handle, "%s\t%s\t%.3f\n",
            proteins->items[edges->items[i].seq_a],
            proteins->items[edges->items[i].seq_b],
            round_score(edges->items[i].normalized_score));
  }
  fclose(handle);
}

static void write_inparalogs_raw(const char *path, EdgeVec *edges, StringVec *proteins, StringVec *taxa) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open inparalog raw output");
  for (size_t i = 0; i < edges->len; i++) {
    fprintf(handle, "%s\t%s\t%s\t%.12f\n",
            proteins->items[edges->items[i].seq_a],
            proteins->items[edges->items[i].seq_b],
            taxa->items[edges->items[i].taxon_id],
            edges->items[i].unnormalized_score);
  }
  fclose(handle);
}

static void write_summary(FILE *handle, size_t record_count, size_t query_group_count, size_t edge_count, size_t protein_count, size_t taxon_count) {
  fprintf(handle, "records\t%zu\n", record_count);
  fprintf(handle, "query_groups\t%zu\n", query_group_count);
  fprintf(handle, "inparalog_edges\t%zu\n", edge_count);
  fprintf(handle, "proteins\t%zu\n", protein_count);
  fprintf(handle, "taxa\t%zu\n", taxon_count);
}

static void free_strings(StringVec *values) {
  for (size_t i = 0; i < values->len; i++) free(values->items[i]);
  free(values->items);
}

static void usage(void) {
  fprintf(stderr, "usage: indexed_inparalogs similarities.bin proteins.tsv taxa.tsv orthologs.txt out_dir percent_match_cutoff evalue_exp_cutoff [shard_index shard_count]\n");
}

int main(int argc, char **argv) {
  if (argc != 8 && argc != 10) {
    usage();
    return 1;
  }

  const char *binary_path = argv[1];
  const char *proteins_path = argv[2];
  const char *taxa_path = argv[3];
  const char *orthologs_path = argv[4];
  const char *out_dir = argv[5];
  float percent_match_cutoff = (float)atof(argv[6]);
  int32_t evalue_exp_cutoff = (int32_t)atoi(argv[7]);
  uint32_t shard_index = 0;
  uint32_t shard_count = 1;
  if (argc == 10) {
    shard_index = (uint32_t)atoi(argv[8]);
    shard_count = (uint32_t)atoi(argv[9]);
    if (shard_count == 0 || shard_index >= shard_count) die("Invalid shard arguments");
  }

  g_records = read_binary_records(binary_path);
  StringVec proteins = read_index_values(proteins_path);
  StringVec taxa = read_index_values(taxa_path);
  RefVec by_query_evalue = build_sorted_refs(g_records.len, cmp_ref_query_evalue);
  RefVec by_query_subject = build_sorted_refs(g_records.len, cmp_ref_query_subject);
  GroupVec groups = build_query_groups(&by_query_evalue);
  NameIndexVec protein_map = build_name_index(&proteins);
  bool *ortholog_mask = load_ortholog_id_mask(orthologs_path, &protein_map, proteins.len);
  EdgeVec edges = build_inparalog_edges(&groups, &by_query_subject, &proteins, ortholog_mask, percent_match_cutoff, evalue_exp_cutoff, shard_index, shard_count);

  make_dir(out_dir);
  char out_path[4096];
  char raw_path[4096];
  char summary_path[4096];
  snprintf(out_path, sizeof(out_path), "%s/inparalogs.indexed.txt", out_dir);
  snprintf(raw_path, sizeof(raw_path), "%s/inparalogs.indexed.raw.tsv", out_dir);
  snprintf(summary_path, sizeof(summary_path), "%s/inparalogs.indexed.summary.tsv", out_dir);
  write_inparalogs(out_path, &edges, &proteins);
  write_inparalogs_raw(raw_path, &edges, &proteins, &taxa);
  FILE *summary = fopen(summary_path, "w");
  if (!summary) die("Could not open summary output");
  write_summary(summary, g_records.len, groups.len, edges.len, proteins.len, taxa.len);
  fclose(summary);

  free(edges.items);
  free(ortholog_mask);
  free(protein_map.items);
  free(groups.items);
  free(by_query_evalue.items);
  free(by_query_subject.items);
  free(g_records.items);
  free_strings(&proteins);
  free_strings(&taxa);
  return 0;
}
