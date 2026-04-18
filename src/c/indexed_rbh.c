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
  char **items;
  size_t len;
  size_t cap;
} StringVec;

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
  uint32_t query_id;
  uint32_t query_taxon_id;
  uint32_t subject_taxon_id;
  uint32_t best_subject_id;
  int32_t best_evalue_exp;
  float best_evalue_mant;
  bool unique_best;
} BestHitGroup;

typedef struct {
  BestHitGroup *items;
  size_t len;
  size_t cap;
} BestHitVec;

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

static void push_best_hit(BestHitVec *vec, BestHitGroup item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 256;
    vec->items = (BestHitGroup *)xrealloc(vec->items, vec->cap * sizeof(BestHitGroup));
  }
  vec->items[vec->len++] = item;
}

static int cmp_best_hit_key(const void *left, const void *right) {
  const BestHitGroup *a = (const BestHitGroup *)left;
  const BestHitGroup *b = (const BestHitGroup *)right;
  if (a->query_id != b->query_id) return a->query_id < b->query_id ? -1 : 1;
  if (a->subject_taxon_id != b->subject_taxon_id) return a->subject_taxon_id < b->subject_taxon_id ? -1 : 1;
  return 0;
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

static RefVec build_sorted_refs(size_t record_count) {
  RefVec refs = {0};
  refs.items = (RecordRef *)xmalloc(record_count * sizeof(RecordRef));
  refs.len = record_count;
  refs.cap = record_count;
  for (size_t i = 0; i < record_count; i++) refs.items[i].record_index = (uint32_t)i;
  qsort(refs.items, refs.len, sizeof(RecordRef), cmp_ref_query_taxon_best);
  return refs;
}

static bool passes_evalue_cutoff(const SimilarityRecordBin *record, float cutoff_mant, int32_t cutoff_exp) {
  if (record->evalue_mant == 0.0f) return true;
  if (record->evalue_exp != cutoff_exp) return record->evalue_exp < cutoff_exp;
  return record->evalue_mant <= cutoff_mant;
}

static BestHitVec build_best_hits(RefVec *sorted_refs, float percent_match_cutoff, float cutoff_mant, int32_t cutoff_exp) {
  BestHitVec groups = {0};
  size_t i = 0;
  while (i < sorted_refs->len) {
    if (i > 0 && i % 1000000 == 0) {
      fprintf(stderr, "[indexed_rbh] processed %zu/%zu sorted similarity references, %zu best-hit buckets so far\n",
              i, sorted_refs->len, groups.len);
    }
    const SimilarityRecordBin *first = &g_records.items[sorted_refs->items[i].record_index];
    uint32_t query_id = first->query_id;
    uint32_t query_taxon_id = first->query_taxon_id;
    uint32_t subject_taxon_id = first->subject_taxon_id;
    size_t end = i;
    BestHitGroup group;
    group.query_id = query_id;
    group.query_taxon_id = query_taxon_id;
    group.subject_taxon_id = subject_taxon_id;
    group.best_subject_id = 0;
    group.best_evalue_exp = 0;
    group.best_evalue_mant = 0.0f;
    group.unique_best = false;
    bool have_best = false;

    while (end < sorted_refs->len) {
      const SimilarityRecordBin *record = &g_records.items[sorted_refs->items[end].record_index];
      if (record->query_id != query_id || record->subject_taxon_id != subject_taxon_id) break;
      if (record->query_taxon_id != record->subject_taxon_id &&
          record->percent_match >= percent_match_cutoff &&
          passes_evalue_cutoff(record, cutoff_mant, cutoff_exp)) {
        if (!have_best) {
          group.best_subject_id = record->subject_id;
          group.best_evalue_exp = record->evalue_exp;
          group.best_evalue_mant = record->evalue_mant;
          group.unique_best = true;
          have_best = true;
        } else if (record->evalue_exp == group.best_evalue_exp &&
                   record->evalue_mant == group.best_evalue_mant &&
                   record->subject_id != group.best_subject_id) {
          group.unique_best = false;
        }
      }
      end++;
    }

    if (have_best) push_best_hit(&groups, group);
    i = end;
  }
  qsort(groups.items, groups.len, sizeof(BestHitGroup), cmp_best_hit_key);
  return groups;
}

static BestHitGroup *find_best_hit(BestHitVec *vec, uint32_t query_id, uint32_t subject_taxon_id) {
  size_t left = 0;
  size_t right = vec->len;
  while (left < right) {
    size_t mid = left + (right - left) / 2;
    BestHitGroup *group = &vec->items[mid];
    if (group->query_id < query_id ||
        (group->query_id == query_id && group->subject_taxon_id < subject_taxon_id)) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  if (left < vec->len &&
      vec->items[left].query_id == query_id &&
      vec->items[left].subject_taxon_id == subject_taxon_id) {
    return &vec->items[left];
  }
  return NULL;
}

static void make_dir(const char *path) {
  char command[4096];
  snprintf(command, sizeof(command), "mkdir -p \"%s\"", path);
  if (system(command) != 0) die("Could not create output directory");
}

static void write_rbh(const char *path, BestHitVec *best_hits, StringVec *proteins, uint32_t shard_index, uint32_t shard_count) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open RBH output");
  for (size_t i = 0; i < best_hits->len; i++) {
    if (i > 0 && i % 250000 == 0) {
      fprintf(stderr, "[indexed_rbh] shard %u/%u wrote %zu/%zu best-hit buckets\n",
              shard_index + 1, shard_count, i, best_hits->len);
    }
    BestHitGroup *left = &best_hits->items[i];
    if (shard_count > 1 && (left->query_id % shard_count) != shard_index) continue;
    if (!left->unique_best) continue;
    BestHitGroup *right = find_best_hit(best_hits, left->best_subject_id, left->query_taxon_id);
    if (!right || !right->unique_best) continue;
    if (right->best_subject_id != left->query_id) continue;
    if (strcmp(proteins->items[left->query_id], proteins->items[left->best_subject_id]) >= 0) continue;
    fprintf(handle, "%s\t%s\n", proteins->items[left->query_id], proteins->items[left->best_subject_id]);
  }
  fclose(handle);
}

static void write_summary(FILE *handle, size_t record_count, size_t best_hit_count, size_t protein_count, size_t taxon_count) {
  fprintf(handle, "records\t%zu\n", record_count);
  fprintf(handle, "best_hit_buckets\t%zu\n", best_hit_count);
  fprintf(handle, "proteins\t%zu\n", protein_count);
  fprintf(handle, "taxa\t%zu\n", taxon_count);
}

static void free_strings(StringVec *values) {
  for (size_t i = 0; i < values->len; i++) free(values->items[i]);
  free(values->items);
}

static void usage(void) {
  fprintf(stderr, "usage: indexed_rbh similarities.bin proteins.tsv taxa.tsv out_dir percent_match_cutoff evalue_cutoff_mant evalue_cutoff_exp [shard_index shard_count]\n");
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
  RefVec refs = build_sorted_refs(g_records.len);
  BestHitVec best_hits = build_best_hits(&refs, percent_match_cutoff, cutoff_mant, cutoff_exp);

  make_dir(out_dir);
  char out_path[4096];
  char summary_path[4096];
  snprintf(out_path, sizeof(out_path), "%s/rbh_1to1.txt", out_dir);
  snprintf(summary_path, sizeof(summary_path), "%s/rbh_1to1.summary.tsv", out_dir);
  write_rbh(out_path, &best_hits, &proteins, shard_index, shard_count);
  FILE *summary = fopen(summary_path, "w");
  if (!summary) die("Could not open summary output");
  write_summary(summary, g_records.len, best_hits.len, proteins.len, taxa.len);
  fclose(summary);

  free(best_hits.items);
  free(refs.items);
  free(g_records.items);
  free_strings(&proteins);
  free_strings(&taxa);
  return 0;
}
