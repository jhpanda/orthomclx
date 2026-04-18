#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

typedef struct {
  char *key;
  uint32_t value;
  int in_use;
} HashEntry;

typedef struct {
  HashEntry *entries;
  size_t capacity;
  size_t size;
} HashMap;

typedef struct {
  char *gene_id;
  char *taxon;
  int length;
} GeneRecord;

typedef struct {
  GeneRecord *records;
  size_t capacity;
  size_t size;
  HashMap by_gene;
} GeneTable;

typedef struct {
  char **items;
  size_t size;
  size_t capacity;
  HashMap index;
} StringIndex;

typedef struct {
  int start;
  int end;
} Span;

typedef struct {
  Span *items;
  size_t size;
  size_t capacity;
} SpanVec;

typedef struct {
  char *query_id;
  char *subject_id;
  uint32_t query_idx;
  uint32_t subject_idx;
  uint32_t query_taxon_idx;
  uint32_t subject_taxon_idx;
  int query_length;
  int subject_length;
  int query_shorter;
  float evalue_mant;
  int32_t evalue_exp;
  double total_identities;
  int total_length;
  SpanVec spans;
  int active;
} SubjectAccumulator;

#pragma pack(push, 1)
typedef struct {
  uint32_t query_idx;
  uint32_t subject_idx;
  uint32_t query_taxon_idx;
  uint32_t subject_taxon_idx;
  float evalue_mant;
  int32_t evalue_exp;
  float percent_identity;
  float percent_match;
} SimilarityRecordPacked;
#pragma pack(pop)

static void die(const char *message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
}

static void report_parser_progress(uint64_t line_number, uint64_t record_count, size_t protein_count, size_t taxon_count) {
  fprintf(stderr,
          "[parse_blast_compiled] processed %llu BLAST lines, wrote %llu similarity records, indexed %zu proteins across %zu taxa\n",
          (unsigned long long) line_number,
          (unsigned long long) record_count,
          protein_count,
          taxon_count);
}

static int is_split_space(char ch) {
  return ch == ' ' || ch == '\t';
}

static int split_whitespace_fields(char *line, char **fields, int max_fields) {
  int count = 0;
  char *cursor = line;
  while (*cursor) {
    while (*cursor && is_split_space(*cursor)) {
      cursor++;
    }
    if (!*cursor) {
      break;
    }
    if (count >= max_fields) {
      return max_fields + 1;
    }
    fields[count++] = cursor;
    while (*cursor && !is_split_space(*cursor)) {
      cursor++;
    }
    if (*cursor) {
      *cursor = '\0';
      cursor++;
    }
  }
  return count;
}

static void *xmalloc(size_t size) {
  void *ptr = malloc(size);
  if (ptr == NULL) {
    die("Out of memory");
  }
  return ptr;
}

static void *xcalloc(size_t count, size_t size) {
  void *ptr = calloc(count, size);
  if (ptr == NULL) {
    die("Out of memory");
  }
  return ptr;
}

static void *xrealloc(void *ptr, size_t size) {
  void *result = realloc(ptr, size);
  if (result == NULL) {
    die("Out of memory");
  }
  return result;
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
      if (length == 0) {
        return 0;
      }
      break;
    }
    if (length + 1 >= *capacity) {
      *capacity *= 2;
      *buffer = (char *)xrealloc(*buffer, *capacity);
    }
    (*buffer)[length++] = (char)ch;
    if (ch == '\n') {
      break;
    }
  }
  (*buffer)[length] = '\0';
  return 1;
}

static char *xstrdup(const char *value) {
  size_t length = strlen(value);
  char *copy = (char *) xmalloc(length + 1);
  memcpy(copy, value, length + 1);
  return copy;
}

static uint64_t hash_string(const char *value) {
  uint64_t hash = 1469598103934665603ULL;
  while (*value) {
    hash ^= (unsigned char) *value++;
    hash *= 1099511628211ULL;
  }
  return hash;
}

static void hashmap_init(HashMap *map, size_t capacity) {
  map->capacity = capacity;
  map->size = 0;
  map->entries = (HashEntry *) xcalloc(capacity, sizeof(HashEntry));
}

static void hashmap_free(HashMap *map) {
  free(map->entries);
  map->entries = NULL;
  map->capacity = 0;
  map->size = 0;
}

static void hashmap_grow(HashMap *map) {
  HashEntry *old_entries = map->entries;
  size_t old_capacity = map->capacity;
  hashmap_init(map, old_capacity * 2);
  for (size_t i = 0; i < old_capacity; i++) {
    if (!old_entries[i].in_use) {
      continue;
    }
    uint64_t hash = hash_string(old_entries[i].key);
    size_t index = (size_t) (hash % map->capacity);
    while (map->entries[index].in_use) {
      index = (index + 1) % map->capacity;
    }
    map->entries[index] = old_entries[i];
    map->size++;
  }
  free(old_entries);
}

static int hashmap_get(const HashMap *map, const char *key, uint32_t *value) {
  if (map->capacity == 0) {
    return 0;
  }
  uint64_t hash = hash_string(key);
  size_t index = (size_t) (hash % map->capacity);
  for (;;) {
    HashEntry *entry = &map->entries[index];
    if (!entry->in_use) {
      return 0;
    }
    if (strcmp(entry->key, key) == 0) {
      if (value) {
        *value = entry->value;
      }
      return 1;
    }
    index = (index + 1) % map->capacity;
  }
}

static void hashmap_put(HashMap *map, char *key_owned, uint32_t value) {
  if ((map->size + 1) * 10 >= map->capacity * 7) {
    hashmap_grow(map);
  }
  uint64_t hash = hash_string(key_owned);
  size_t index = (size_t) (hash % map->capacity);
  while (map->entries[index].in_use) {
    if (strcmp(map->entries[index].key, key_owned) == 0) {
      map->entries[index].value = value;
      free(key_owned);
      return;
    }
    index = (index + 1) % map->capacity;
  }
  map->entries[index].key = key_owned;
  map->entries[index].value = value;
  map->entries[index].in_use = 1;
  map->size++;
}

static void string_index_init(StringIndex *index) {
  index->items = NULL;
  index->size = 0;
  index->capacity = 0;
  hashmap_init(&index->index, 1024);
}

static void string_index_free(StringIndex *index) {
  for (size_t i = 0; i < index->size; i++) {
    free(index->items[i]);
  }
  free(index->items);
  hashmap_free(&index->index);
}

static uint32_t string_index_get_or_add(StringIndex *index, const char *value) {
  uint32_t existing = 0;
  if (hashmap_get(&index->index, value, &existing)) {
    return existing;
  }
  if (index->size == index->capacity) {
    index->capacity = index->capacity == 0 ? 64 : index->capacity * 2;
    index->items = (char **) xrealloc(index->items, index->capacity * sizeof(char *));
  }
  uint32_t new_index = (uint32_t) index->size;
  char *copy = xstrdup(value);
  index->items[index->size++] = copy;
  hashmap_put(&index->index, xstrdup(value), new_index);
  return new_index;
}

static void gene_table_init(GeneTable *table) {
  table->records = NULL;
  table->size = 0;
  table->capacity = 0;
  hashmap_init(&table->by_gene, 2048);
}

static void gene_table_free(GeneTable *table) {
  for (size_t i = 0; i < table->size; i++) {
    free(table->records[i].gene_id);
    free(table->records[i].taxon);
  }
  free(table->records);
  hashmap_free(&table->by_gene);
}

static void gene_table_add(GeneTable *table, const char *gene_id, const char *taxon, int length) {
  uint32_t existing = 0;
  if (hashmap_get(&table->by_gene, gene_id, &existing)) {
    die("Duplicate gene ID found in FASTA input");
  }
  if (table->size == table->capacity) {
    table->capacity = table->capacity == 0 ? 1024 : table->capacity * 2;
    table->records = (GeneRecord *) xrealloc(table->records, table->capacity * sizeof(GeneRecord));
  }
  uint32_t index = (uint32_t) table->size;
  table->records[index].gene_id = xstrdup(gene_id);
  table->records[index].taxon = xstrdup(taxon);
  table->records[index].length = length;
  table->size++;
  hashmap_put(&table->by_gene, xstrdup(gene_id), index);
}

static const GeneRecord *gene_table_get(const GeneTable *table, const char *gene_id) {
  uint32_t index = 0;
  if (!hashmap_get(&table->by_gene, gene_id, &index)) {
    return NULL;
  }
  return &table->records[index];
}

static void span_vec_init(SpanVec *vec) {
  vec->items = NULL;
  vec->size = 0;
  vec->capacity = 0;
}

static void span_vec_free(SpanVec *vec) {
  free(vec->items);
  vec->items = NULL;
  vec->size = 0;
  vec->capacity = 0;
}

static void span_vec_clear(SpanVec *vec) {
  vec->size = 0;
}

static void span_vec_push(SpanVec *vec, int start, int end) {
  if (vec->size == vec->capacity) {
    vec->capacity = vec->capacity == 0 ? 16 : vec->capacity * 2;
    vec->items = (Span *) xrealloc(vec->items, vec->capacity * sizeof(Span));
  }
  vec->items[vec->size].start = start;
  vec->items[vec->size].end = end;
  vec->size++;
}

static int span_compare(const void *left, const void *right) {
  const Span *a = (const Span *) left;
  const Span *b = (const Span *) right;
  int a_start = a->start < a->end ? a->start : a->end;
  int b_start = b->start < b->end ? b->start : b->end;
  if (a_start < b_start) return -1;
  if (a_start > b_start) return 1;
  return 0;
}

static int compute_non_overlapping_match_length(SpanVec *vec) {
  if (vec->size == 0) {
    return 0;
  }
  qsort(vec->items, vec->size, sizeof(Span), span_compare);
  int start = vec->items[0].start < vec->items[0].end ? vec->items[0].start : vec->items[0].end;
  int end = vec->items[0].start < vec->items[0].end ? vec->items[0].end : vec->items[0].start;
  int length = 0;
  for (size_t i = 1; i < vec->size; i++) {
    int hsp_start = vec->items[i].start < vec->items[i].end ? vec->items[i].start : vec->items[i].end;
    int hsp_end = vec->items[i].start < vec->items[i].end ? vec->items[i].end : vec->items[i].start;
    if (hsp_end <= end) {
      continue;
    }
    if (hsp_start <= end) {
      end = hsp_end;
    } else {
      length += end - start + 1;
      start = hsp_start;
      end = hsp_end;
    }
  }
  length += end - start + 1;
  return length;
}

static void ensure_parent_dir(const char *path) {
  char buffer[4096];
  size_t length = strlen(path);
  if (length >= sizeof(buffer)) {
    die("Path too long");
  }
  memcpy(buffer, path, length + 1);
  for (size_t i = 1; i < length; i++) {
    if (buffer[i] == '/') {
      buffer[i] = '\0';
      mkdir(buffer, 0777);
      buffer[i] = '/';
    }
  }
  mkdir(buffer, 0777);
}

static char *trim_line(char *line) {
  size_t length = strlen(line);
  while (length > 0 && (line[length - 1] == '\n' || line[length - 1] == '\r')) {
    line[--length] = '\0';
  }
  return line;
}

static void read_fasta_lengths(const char *fasta_dir, GeneTable *genes) {
  DIR *directory = opendir(fasta_dir);
  if (directory == NULL) {
    die("Could not open FASTA directory");
  }
  struct dirent *entry = NULL;
  while ((entry = readdir(directory)) != NULL) {
    const char *name = entry->d_name;
    size_t name_len = strlen(name);
    if (name[0] == '.') {
      continue;
    }
    if (name_len < 7 || strcmp(name + name_len - 6, ".fasta") != 0) {
      continue;
    }
    char path[4096];
    snprintf(path, sizeof(path), "%s/%s", fasta_dir, name);
    FILE *handle = fopen(path, "r");
    if (handle == NULL) {
      closedir(directory);
      die("Could not open FASTA file");
    }
    char taxon[256];
    if (name_len - 6 >= sizeof(taxon)) {
      fclose(handle);
      closedir(directory);
      die("Taxon name too long");
    }
    memcpy(taxon, name, name_len - 6);
    taxon[name_len - 6] = '\0';

    char *line = NULL;
    size_t line_capacity = 0;
    char *gene_id = NULL;
    int gene_length = 0;
    while (read_line_alloc(handle, &line, &line_capacity)) {
      trim_line(line);
      if (line[0] == '\0') {
        continue;
      }
      if (line[0] == '>') {
        if (gene_id != NULL) {
          gene_table_add(genes, gene_id, taxon, gene_length);
          free(gene_id);
        }
        char *cursor = line + 1;
        while (*cursor && isspace((unsigned char) *cursor)) cursor++;
        char *end = cursor;
        while (*end && !isspace((unsigned char) *end)) end++;
        char saved = *end;
        *end = '\0';
        gene_id = xstrdup(cursor);
        *end = saved;
        gene_length = 0;
      } else {
        gene_length += (int) strlen(line);
      }
    }
    if (gene_id != NULL) {
      gene_table_add(genes, gene_id, taxon, gene_length);
      free(gene_id);
    }
    free(line);
    fclose(handle);
  }
  closedir(directory);
}

static void format_evalue(const char *evalue_text, float *mant_out, int32_t *exp_out) {
  char normalized[128];
  if (evalue_text[0] == 'e' || evalue_text[0] == 'E') {
    snprintf(normalized, sizeof(normalized), "1%s", evalue_text);
  } else {
    snprintf(normalized, sizeof(normalized), "%s", evalue_text);
  }
  double value = strtod(normalized, NULL);
  char scientific[128];
  snprintf(scientific, sizeof(scientific), "%.3e", value);
  char *e_ptr = strchr(scientific, 'e');
  if (e_ptr == NULL) {
    die("Could not parse e-value");
  }
  *e_ptr = '\0';
  double mantissa = strtod(scientific, NULL);
  int exponent = (int) strtol(e_ptr + 1, NULL, 10);
  mantissa = floor(mantissa * 100.0 + 0.5) / 100.0;
  *mant_out = (float) mantissa;
  *exp_out = exponent;
}

static float round_tenth(double value) {
  return (float) (floor(value * 10.0 + 0.5) / 10.0);
}

static void accumulator_init(SubjectAccumulator *acc) {
  memset(acc, 0, sizeof(*acc));
  span_vec_init(&acc->spans);
}

static void accumulator_reset(SubjectAccumulator *acc) {
  free(acc->query_id);
  free(acc->subject_id);
  acc->query_id = NULL;
  acc->subject_id = NULL;
  acc->active = 0;
  span_vec_clear(&acc->spans);
  acc->total_identities = 0.0;
  acc->total_length = 0;
}

static void accumulator_free(SubjectAccumulator *acc) {
  accumulator_reset(acc);
  span_vec_free(&acc->spans);
}

static void accumulator_start(
  SubjectAccumulator *acc,
  const char *query_id,
  const char *subject_id,
  uint32_t query_idx,
  uint32_t subject_idx,
  uint32_t query_taxon_idx,
  uint32_t subject_taxon_idx,
  const GeneRecord *query_gene,
  const GeneRecord *subject_gene,
  const char *evalue_text
) {
  accumulator_reset(acc);
  acc->query_id = xstrdup(query_id);
  acc->subject_id = xstrdup(subject_id);
  acc->query_idx = query_idx;
  acc->subject_idx = subject_idx;
  acc->query_taxon_idx = query_taxon_idx;
  acc->subject_taxon_idx = subject_taxon_idx;
  acc->query_length = query_gene->length;
  acc->subject_length = subject_gene->length;
  acc->query_shorter = query_gene->length < subject_gene->length;
  format_evalue(evalue_text, &acc->evalue_mant, &acc->evalue_exp);
  acc->active = 1;
}

static void write_accumulator_record(
  SubjectAccumulator *acc,
  FILE *binary_handle,
  uint64_t *record_count,
  int32_t *min_nonzero_exp
) {
  if (!acc->active) {
    return;
  }
  int non_overlap = compute_non_overlapping_match_length(&acc->spans);
  int shorter_length = acc->query_shorter ? acc->query_length : acc->subject_length;
  SimilarityRecordPacked record;
  record.query_idx = acc->query_idx;
  record.subject_idx = acc->subject_idx;
  record.query_taxon_idx = acc->query_taxon_idx;
  record.subject_taxon_idx = acc->subject_taxon_idx;
  record.evalue_mant = acc->evalue_mant;
  record.evalue_exp = acc->evalue_exp;
  record.percent_identity = round_tenth(acc->total_identities / (double) acc->total_length);
  record.percent_match = round_tenth(non_overlap * 100.0 / (double) shorter_length);
  if (record.evalue_mant != 0.0f) {
    if (*min_nonzero_exp == 0 || record.evalue_exp < *min_nonzero_exp) {
      *min_nonzero_exp = record.evalue_exp;
    }
  }
  if (fwrite(&record, sizeof(record), 1, binary_handle) != 1) {
    die("Could not write compiled similarity record");
  }
  (*record_count)++;
}

static void rewrite_zero_evalues(const char *binary_path, int32_t adjusted_zero_exp) {
  FILE *handle = fopen(binary_path, "r+b");
  if (handle == NULL) {
    die("Could not reopen binary file for zero e-value rewrite");
  }
  SimilarityRecordPacked record;
  while (fread(&record, sizeof(record), 1, handle) == 1) {
    if (record.evalue_mant == 0.0f && record.evalue_exp == 0) {
      record.evalue_exp = adjusted_zero_exp;
      if (fseek(handle, -(long) sizeof(record), SEEK_CUR) != 0) {
        fclose(handle);
        die("Could not seek while rewriting binary file");
      }
      if (fwrite(&record, sizeof(record), 1, handle) != 1) {
        fclose(handle);
        die("Could not rewrite zero e-value exponent");
      }
    }
  }
  fclose(handle);
}

static void write_index_tsv(const char *path, const StringIndex *index) {
  FILE *handle = fopen(path, "w");
  if (handle == NULL) {
    die("Could not open TSV index file");
  }
  for (size_t i = 0; i < index->size; i++) {
    fprintf(handle, "%zu\t%s\n", i, index->items[i]);
  }
  fclose(handle);
}

static void parse_blast_to_compiled(const char *blast_path, const char *fasta_dir, const char *out_dir) {
  GeneTable genes;
  gene_table_init(&genes);
  read_fasta_lengths(fasta_dir, &genes);

  StringIndex proteins;
  StringIndex taxa;
  string_index_init(&proteins);
  string_index_init(&taxa);

  ensure_parent_dir(out_dir);
  mkdir(out_dir, 0777);

  char binary_path[4096];
  char proteins_path[4096];
  char taxa_path[4096];
  snprintf(binary_path, sizeof(binary_path), "%s/similarities.bin", out_dir);
  snprintf(proteins_path, sizeof(proteins_path), "%s/proteins.tsv", out_dir);
  snprintf(taxa_path, sizeof(taxa_path), "%s/taxa.tsv", out_dir);

  FILE *blast_handle = fopen(blast_path, "r");
  if (blast_handle == NULL) {
    die("Could not open BLAST file");
  }
  FILE *binary_handle = fopen(binary_path, "wb");
  if (binary_handle == NULL) {
    fclose(blast_handle);
    die("Could not open compiled binary output");
  }

  SubjectAccumulator acc;
  accumulator_init(&acc);
  uint64_t record_count = 0;
  int32_t min_nonzero_exp = 0;

  char *line = NULL;
  size_t line_capacity = 0;
  uint64_t line_number = 0;
  while (read_line_alloc(blast_handle, &line, &line_capacity)) {
    line_number++;
    trim_line(line);
    if (line[0] == '\0') {
      continue;
    }

    char *fields[12];
    int field_count = 0;
    field_count = split_whitespace_fields(line, fields, 12);
    if (field_count != 12) {
      fprintf(stderr, "BLAST line %llu does not have 12 columns\n", (unsigned long long) line_number);
      exit(1);
    }

    const char *query_id = fields[0];
    const char *subject_id = fields[1];
    const char *percent_identity = fields[2];
    const char *length_text = fields[3];
    const char *query_start = fields[6];
    const char *query_end = fields[7];
    const char *subject_start = fields[8];
    const char *subject_end = fields[9];
    const char *evalue_text = fields[10];

    const GeneRecord *query_gene = gene_table_get(&genes, query_id);
    const GeneRecord *subject_gene = gene_table_get(&genes, subject_id);
    if (query_gene == NULL || subject_gene == NULL) {
      die("BLAST referenced a protein not found in FASTA input");
    }

    uint32_t query_idx = string_index_get_or_add(&proteins, query_id);
    uint32_t subject_idx = string_index_get_or_add(&proteins, subject_id);
    uint32_t query_taxon_idx = string_index_get_or_add(&taxa, query_gene->taxon);
    uint32_t subject_taxon_idx = string_index_get_or_add(&taxa, subject_gene->taxon);

    if (!acc.active || strcmp(query_id, acc.query_id) != 0 || strcmp(subject_id, acc.subject_id) != 0) {
      write_accumulator_record(&acc, binary_handle, &record_count, &min_nonzero_exp);
      accumulator_start(
        &acc,
        query_id,
        subject_id,
        query_idx,
        subject_idx,
        query_taxon_idx,
        subject_taxon_idx,
        query_gene,
        subject_gene,
        evalue_text
      );
    }

    int start = acc.query_shorter ? atoi(query_start) : atoi(subject_start);
    int end = acc.query_shorter ? atoi(query_end) : atoi(subject_end);
    span_vec_push(&acc.spans, start, end);
    acc.total_identities += atof(percent_identity) * atoi(length_text);
    acc.total_length += atoi(length_text);

    if (line_number % 1000000ULL == 0) {
      report_parser_progress(line_number, record_count, proteins.size, taxa.size);
    }
  }

  write_accumulator_record(&acc, binary_handle, &record_count, &min_nonzero_exp);
  free(line);
  accumulator_free(&acc);
  fclose(binary_handle);
  fclose(blast_handle);

  if (min_nonzero_exp != 0) {
    rewrite_zero_evalues(binary_path, min_nonzero_exp - 1);
  }

  write_index_tsv(proteins_path, &proteins);
  write_index_tsv(taxa_path, &taxa);

  size_t protein_count = proteins.size;
  size_t taxon_count = taxa.size;

  string_index_free(&proteins);
  string_index_free(&taxa);
  gene_table_free(&genes);

  fprintf(stdout, "compiled %llu records, %zu proteins, %zu taxa\n",
    (unsigned long long) record_count, protein_count, taxon_count);
}

int main(int argc, char **argv) {
  if (argc != 4) {
    fprintf(stderr, "Usage: %s <blast.tsv> <fasta_dir> <out_dir>\n", argv[0]);
    return 1;
  }
  parse_blast_to_compiled(argv[1], argv[2], argv[3]);
  return 0;
}
