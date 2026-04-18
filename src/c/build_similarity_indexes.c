#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

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
  SimilarityRecordBin *items;
  size_t len;
  size_t cap;
} RecordVec;

typedef struct {
  RecordRef *items;
  size_t len;
  size_t cap;
} RefVec;

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

static void push_record(RecordVec *vec, SimilarityRecordBin item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 1024;
    vec->items = (SimilarityRecordBin *)xrealloc(vec->items, vec->cap * sizeof(SimilarityRecordBin));
  }
  vec->items[vec->len++] = item;
}

static RecordVec read_binary_records(const char *path) {
  RecordVec records = {0};
  FILE *handle = fopen(path, "rb");
  if (!handle) die("Could not open similarities.bin");
  SimilarityRecordBin record;
  while (fread(&record, sizeof(SimilarityRecordBin), 1, handle) == 1) {
    push_record(&records, record);
  }
  fclose(handle);
  return records;
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

static RefVec build_sorted_refs(size_t record_count, int (*cmp)(const void *, const void *), const char *label) {
  RefVec refs = {0};
  refs.len = record_count;
  refs.cap = record_count;
  refs.items = (RecordRef *)xmalloc((record_count ? record_count : 1) * sizeof(RecordRef));
  for (size_t i = 0; i < record_count; i++) {
    refs.items[i].record_index = (uint32_t)i;
  }
  fprintf(stderr, "[build_similarity_indexes] sorting %s refs (%zu records)\n", label, record_count);
  qsort(refs.items, refs.len, sizeof(RecordRef), cmp);
  return refs;
}

static void write_refs(const char *path, RefVec *refs) {
  FILE *handle = fopen(path, "wb");
  if (!handle) die("Could not open output ref index file");
  if (refs->len && fwrite(refs->items, sizeof(RecordRef), refs->len, handle) != refs->len) {
    fclose(handle);
    die("Could not write ref index file");
  }
  fclose(handle);
}

static void join_path(char *out, size_t out_size, const char *dir, const char *name) {
  if (snprintf(out, out_size, "%s/%s", dir, name) >= (int)out_size) {
    die("Output path too long");
  }
}

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: build_similarity_indexes similarities.bin out_dir\n");
    return 1;
  }

  const char *binary_path = argv[1];
  const char *out_dir = argv[2];
  char query_taxon_path[4096];
  char query_subject_path[4096];
  char query_evalue_path[4096];

  g_records = read_binary_records(binary_path);
  fprintf(stderr, "[build_similarity_indexes] loaded %zu similarity records\n", g_records.len);

  RefVec by_query_taxon = build_sorted_refs(g_records.len, cmp_ref_query_taxon_best, "query_taxon_best");
  RefVec by_query_subject = build_sorted_refs(g_records.len, cmp_ref_query_subject, "query_subject");
  RefVec by_query_evalue = build_sorted_refs(g_records.len, cmp_ref_query_evalue, "query_evalue");

  join_path(query_taxon_path, sizeof(query_taxon_path), out_dir, "refs.query_taxon_best.bin");
  join_path(query_subject_path, sizeof(query_subject_path), out_dir, "refs.query_subject.bin");
  join_path(query_evalue_path, sizeof(query_evalue_path), out_dir, "refs.query_evalue.bin");

  write_refs(query_taxon_path, &by_query_taxon);
  write_refs(query_subject_path, &by_query_subject);
  write_refs(query_evalue_path, &by_query_evalue);

  free(by_query_taxon.items);
  free(by_query_subject.items);
  free(by_query_evalue.items);
  free(g_records.items);
  fprintf(stderr, "[build_similarity_indexes] wrote reusable reference indexes\n");
  return 0;
}
