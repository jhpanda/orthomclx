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
  uint32_t subject_taxon_id;
  uint32_t start;
  uint32_t end;
  int32_t best_evalue_exp;
  float best_evalue_mant;
} QueryTaxonGroup;

typedef struct {
  uint32_t seq_a;
  uint32_t seq_b;
  uint32_t taxon_a;
  uint32_t taxon_b;
  double unnormalized_score;
  double normalized_score;
} OrthologEdge;

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
  QueryTaxonGroup *items;
  size_t len;
  size_t cap;
} GroupVec;

typedef struct {
  OrthologEdge *items;
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

static void push_record(RecordVec *vec, SimilarityRecordBin item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 1024;
    vec->items = (SimilarityRecordBin *)xrealloc(vec->items, vec->cap * sizeof(SimilarityRecordBin));
  }
  vec->items[vec->len++] = item;
}

static void push_group(GroupVec *vec, QueryTaxonGroup item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 256;
    vec->items = (QueryTaxonGroup *)xrealloc(vec->items, vec->cap * sizeof(QueryTaxonGroup));
  }
  vec->items[vec->len++] = item;
}

static void push_edge(EdgeVec *vec, OrthologEdge item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 256;
    vec->items = (OrthologEdge *)xrealloc(vec->items, vec->cap * sizeof(OrthologEdge));
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

static int cmp_ref_query_taxon_best(const void *left, const void *right) {
  const SimilarityRecordBin *a = &g_records.items[((const RecordRef *)left)->record_index];
  const SimilarityRecordBin *b = &g_records.items[((const RecordRef *)right)->record_index];
  if (a->query_id != b->query_id) return a->query_id < b->query_id ? -1 : 1;
  if (a->subject_taxon_id != b->subject_taxon_id) return a->subject_taxon_id < b->subject_taxon_id ? -1 : 1;
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
  const OrthologEdge *a = (const OrthologEdge *)left;
  const OrthologEdge *b = (const OrthologEdge *)right;
  if (a->taxon_a != b->taxon_a) return a->taxon_a < b->taxon_a ? -1 : 1;
  if (a->taxon_b != b->taxon_b) return a->taxon_b < b->taxon_b ? -1 : 1;
  if (a->seq_a != b->seq_a) return a->seq_a < b->seq_a ? -1 : 1;
  if (a->seq_b != b->seq_b) return a->seq_b < b->seq_b ? -1 : 1;
  return 0;
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

static RefVec build_sorted_refs(size_t record_count, int (*cmp)(const void *, const void *)) {
  RefVec refs = {0};
  refs.items = (RecordRef *)xmalloc(record_count * sizeof(RecordRef));
  refs.len = record_count;
  refs.cap = record_count;
  for (size_t i = 0; i < record_count; i++) refs.items[i].record_index = (uint32_t)i;
  qsort(refs.items, refs.len, sizeof(RecordRef), cmp);
  return refs;
}

static GroupVec build_query_taxon_groups(RefVec *sorted_refs) {
  GroupVec groups = {0};
  size_t i = 0;
  while (i < sorted_refs->len) {
    const SimilarityRecordBin *first = &g_records.items[sorted_refs->items[i].record_index];
    uint32_t query_id = first->query_id;
    uint32_t subject_taxon_id = first->subject_taxon_id;
    size_t start = i;
    size_t end = i;
    int32_t best_evalue_exp = 0;
    float best_evalue_mant = 0.0f;
    bool have_best = false;
    while (end < sorted_refs->len) {
      const SimilarityRecordBin *record = &g_records.items[sorted_refs->items[end].record_index];
      if (record->query_id != query_id || record->subject_taxon_id != subject_taxon_id) break;
      if (record->query_taxon_id != record->subject_taxon_id && !have_best) {
        best_evalue_exp = record->evalue_exp;
        best_evalue_mant = record->evalue_mant;
        have_best = true;
      }
      end++;
    }
    QueryTaxonGroup group;
    group.query_id = query_id;
    group.subject_taxon_id = subject_taxon_id;
    group.start = (uint32_t)start;
    group.end = (uint32_t)end;
    group.best_evalue_exp = have_best ? best_evalue_exp : INT32_MAX;
    group.best_evalue_mant = have_best ? best_evalue_mant : 0.0f;
    push_group(&groups, group);
    i = end;
  }
  return groups;
}

static QueryTaxonGroup *find_group(GroupVec *groups, uint32_t query_id, uint32_t subject_taxon_id) {
  size_t left = 0;
  size_t right = groups->len;
  while (left < right) {
    size_t mid = left + (right - left) / 2;
    QueryTaxonGroup *group = &groups->items[mid];
    if (group->query_id < query_id || (group->query_id == query_id && group->subject_taxon_id < subject_taxon_id)) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  if (left < groups->len &&
      groups->items[left].query_id == query_id &&
      groups->items[left].subject_taxon_id == subject_taxon_id) {
    return &groups->items[left];
  }
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

static bool qualifies_as_best_hit(const SimilarityRecordBin *record, const QueryTaxonGroup *group) {
  if (group->best_evalue_exp == INT32_MAX) return false;
  if (record->evalue_mant < 0.01f) return true;
  return record->evalue_exp == group->best_evalue_exp && record->evalue_mant == group->best_evalue_mant;
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

static EdgeVec build_ortholog_edges(GroupVec *groups, RefVec *by_query_taxon, RefVec *by_query_subject, StringVec *proteins, float percent_match_cutoff, float cutoff_mant, int32_t cutoff_exp, uint32_t shard_index, uint32_t shard_count) {
  EdgeVec edges = {0};
  for (size_t i = 0; i < groups->len; i++) {
    QueryTaxonGroup *group = &groups->items[i];
    if (shard_count > 1 && (group->query_id % shard_count) != shard_index) continue;
    if (i > 0 && i % 50000 == 0) {
      fprintf(stderr, "[indexed_orthologs] shard %u/%u processed %zu/%zu query-taxon groups, %zu orthologs so far\n",
              shard_index + 1, shard_count, i, groups->len, edges.len);
    }
    if (group->best_evalue_exp == INT32_MAX) continue;
    for (uint32_t pos = group->start; pos < group->end; pos++) {
      SimilarityRecordBin *left = &g_records.items[by_query_taxon->items[pos].record_index];
      if (!passes_thresholds(left, percent_match_cutoff, cutoff_mant, cutoff_exp)) continue;
      if (!qualifies_as_best_hit(left, group)) continue;
      if (strcmp(proteins->items[left->query_id], proteins->items[left->subject_id]) >= 0) continue;

      SimilarityRecordBin *right = find_record_by_query_subject(by_query_subject, left->subject_id, left->query_id);
      if (!right) continue;
      if (right->query_taxon_id == right->subject_taxon_id) continue;
      if (!passes_thresholds(right, percent_match_cutoff, cutoff_mant, cutoff_exp)) continue;

      QueryTaxonGroup *reverse_group = find_group(groups, right->query_id, right->subject_taxon_id);
      if (!reverse_group) continue;
      if (!qualifies_as_best_hit(right, reverse_group)) continue;

      OrthologEdge edge;
      if (strcmp(proteins->items[left->query_id], proteins->items[left->subject_id]) < 0) {
        edge.seq_a = left->query_id;
        edge.seq_b = left->subject_id;
        edge.taxon_a = left->query_taxon_id;
        edge.taxon_b = left->subject_taxon_id;
      } else {
        edge.seq_a = left->subject_id;
        edge.seq_b = left->query_id;
        edge.taxon_a = left->subject_taxon_id;
        edge.taxon_b = left->query_taxon_id;
      }
      edge.unnormalized_score = score_from_records(left, right, 0.01);
      edge.normalized_score = 0.0;
      push_edge(&edges, edge);
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

  qsort(edges.items, edges.len, sizeof(OrthologEdge), cmp_edge_output);
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

static void write_orthologs(const char *path, EdgeVec *edges, StringVec *proteins) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open ortholog output");
  for (size_t i = 0; i < edges->len; i++) {
    fprintf(handle, "%s\t%s\t%.3f\n",
            proteins->items[edges->items[i].seq_a],
            proteins->items[edges->items[i].seq_b],
            round_score(edges->items[i].normalized_score));
  }
  fclose(handle);
}

static void write_orthologs_raw(const char *path, EdgeVec *edges, StringVec *proteins, StringVec *taxa) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open ortholog raw output");
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

static void write_summary(FILE *handle, size_t record_count, size_t group_count, size_t edge_count, size_t protein_count, size_t taxon_count) {
  fprintf(handle, "records\t%zu\n", record_count);
  fprintf(handle, "query_taxon_groups\t%zu\n", group_count);
  fprintf(handle, "ortholog_edges\t%zu\n", edge_count);
  fprintf(handle, "proteins\t%zu\n", protein_count);
  fprintf(handle, "taxa\t%zu\n", taxon_count);
}

static void free_strings(StringVec *values) {
  for (size_t i = 0; i < values->len; i++) free(values->items[i]);
  free(values->items);
}

static void usage(void) {
  fprintf(stderr, "usage: indexed_orthologs similarities.bin proteins.tsv taxa.tsv out_dir percent_match_cutoff evalue_cutoff_mant evalue_cutoff_exp [shard_index shard_count]\n");
}

int main(int argc, char **argv) {
  if (argc != 8 && argc != 10) {
    usage();
    return 1;
  }

  const char *binary_path = argv[1];
  const char *proteins_path = argv[2];
  const char *taxa_path = argv[3];
  const char *out_dir = argv[4];
  float percent_match_cutoff = (float)atof(argv[5]);
  float cutoff_mant = (float)atof(argv[6]);
  int32_t cutoff_exp = (int32_t)atoi(argv[7]);
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
  RefVec by_query_taxon = build_sorted_refs(g_records.len, cmp_ref_query_taxon_best);
  RefVec by_query_subject = build_sorted_refs(g_records.len, cmp_ref_query_subject);
  GroupVec groups = build_query_taxon_groups(&by_query_taxon);
  EdgeVec orthologs = build_ortholog_edges(&groups, &by_query_taxon, &by_query_subject, &proteins, percent_match_cutoff, cutoff_mant, cutoff_exp, shard_index, shard_count);

  make_dir(out_dir);
  char orthologs_path[4096];
  char raw_path[4096];
  char summary_path[4096];
  snprintf(orthologs_path, sizeof(orthologs_path), "%s/orthologs.indexed.txt", out_dir);
  snprintf(raw_path, sizeof(raw_path), "%s/orthologs.indexed.raw.tsv", out_dir);
  snprintf(summary_path, sizeof(summary_path), "%s/orthologs.indexed.summary.tsv", out_dir);
  write_orthologs(orthologs_path, &orthologs, &proteins);
  write_orthologs_raw(raw_path, &orthologs, &proteins, &taxa);
  FILE *summary = fopen(summary_path, "w");
  if (!summary) die("Could not open summary output");
  write_summary(summary, g_records.len, groups.len, orthologs.len, proteins.len, taxa.len);
  fclose(summary);

  free(orthologs.items);
  free(groups.items);
  free(by_query_taxon.items);
  free(by_query_subject.items);
  free(g_records.items);
  free_strings(&proteins);
  free_strings(&taxa);
  return 0;
}
