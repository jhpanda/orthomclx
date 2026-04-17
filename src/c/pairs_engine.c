#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char *query_id;
  char *subject_id;
  char *query_taxon;
  char *subject_taxon;
  double evalue_mant;
  int evalue_exp;
  double percent_identity;
  double percent_match;
} SimilarityRecord;

typedef struct {
  char *seq_a;
  char *seq_b;
  char *taxon_a;
  char *taxon_b;
  double unnormalized_score;
  double normalized_score;
} EdgeRecord;

typedef struct {
  SimilarityRecord *items;
  size_t len;
  size_t cap;
} RecordVec;

typedef struct {
  EdgeRecord *items;
  size_t len;
  size_t cap;
} EdgeVec;

typedef struct {
  double *items;
  size_t len;
  size_t cap;
} DoubleVec;

typedef struct {
  char **items;
  size_t len;
  size_t cap;
} StringVec;

static void die(const char *message) {
  fprintf(stderr, "%s\n", message);
  exit(1);
}

static void *xmalloc(size_t size) {
  void *ptr = malloc(size);
  if (!ptr) {
    die("Out of memory");
  }
  return ptr;
}

static void *xrealloc(void *ptr, size_t size) {
  void *out = realloc(ptr, size);
  if (!out) {
    die("Out of memory");
  }
  return out;
}

static char *xstrdup(const char *src) {
  size_t len = strlen(src);
  char *dst = (char *)xmalloc(len + 1);
  memcpy(dst, src, len + 1);
  return dst;
}

static void push_record(RecordVec *vec, SimilarityRecord item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 64;
    vec->items = (SimilarityRecord *)xrealloc(vec->items, vec->cap * sizeof(SimilarityRecord));
  }
  vec->items[vec->len++] = item;
}

static void push_edge(EdgeVec *vec, EdgeRecord item) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 32;
    vec->items = (EdgeRecord *)xrealloc(vec->items, vec->cap * sizeof(EdgeRecord));
  }
  vec->items[vec->len++] = item;
}

static void push_double(DoubleVec *vec, double value) {
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 16;
    vec->items = (double *)xrealloc(vec->items, vec->cap * sizeof(double));
  }
  vec->items[vec->len++] = value;
}

static int string_index(StringVec *vec, const char *value) {
  for (size_t i = 0; i < vec->len; i++) {
    if (strcmp(vec->items[i], value) == 0) {
      return (int)i;
    }
  }
  if (vec->len == vec->cap) {
    vec->cap = vec->cap ? vec->cap * 2 : 16;
    vec->items = (char **)xrealloc(vec->items, vec->cap * sizeof(char *));
  }
  vec->items[vec->len] = xstrdup(value);
  vec->len += 1;
  return (int)(vec->len - 1);
}

static int cmp_record_pair(const void *left, const void *right) {
  const SimilarityRecord *a = (const SimilarityRecord *)left;
  const SimilarityRecord *b = (const SimilarityRecord *)right;
  int cmp = strcmp(a->query_id, b->query_id);
  if (cmp != 0) return cmp;
  return strcmp(a->subject_id, b->subject_id);
}

static int cmp_edge_output(const void *left, const void *right) {
  const EdgeRecord *a = (const EdgeRecord *)left;
  const EdgeRecord *b = (const EdgeRecord *)right;
  int cmp = strcmp(a->taxon_a, b->taxon_a);
  if (cmp != 0) return cmp;
  cmp = strcmp(a->taxon_b, b->taxon_b);
  if (cmp != 0) return cmp;
  cmp = strcmp(a->seq_a, b->seq_a);
  if (cmp != 0) return cmp;
  return strcmp(a->seq_b, b->seq_b);
}

static int cmp_edge_union(const void *left, const void *right) {
  const EdgeRecord *a = (const EdgeRecord *)left;
  const EdgeRecord *b = (const EdgeRecord *)right;
  int cmp = strcmp(a->seq_a, b->seq_a);
  if (cmp != 0) return cmp;
  return strcmp(a->seq_b, b->seq_b);
}

static int compare_evalue(const SimilarityRecord *a, const SimilarityRecord *b) {
  if (a->evalue_exp != b->evalue_exp) {
    return a->evalue_exp < b->evalue_exp ? -1 : 1;
  }
  if (a->evalue_mant < b->evalue_mant) return -1;
  if (a->evalue_mant > b->evalue_mant) return 1;
  return 0;
}

static bool passes_thresholds(const SimilarityRecord *record, double percent_match_cutoff, int evalue_exp_cutoff) {
  return record->percent_match >= percent_match_cutoff && record->evalue_exp <= evalue_exp_cutoff;
}

static double score_from_records(const SimilarityRecord *left, const SimilarityRecord *right, double zero_cutoff) {
  if (left->evalue_mant < zero_cutoff || right->evalue_mant < zero_cutoff) {
    return (left->evalue_exp + right->evalue_exp) / -2.0;
  }
  return (log10(left->evalue_mant * right->evalue_mant) + left->evalue_exp + right->evalue_exp) / -2.0;
}

static double rounded_score(double score) {
  return floor(score * 1000.0 + 0.5) / 1000.0;
}

static SimilarityRecord *find_pair(RecordVec *records, const char *query_id, const char *subject_id) {
  size_t left = 0;
  size_t right = records->len;
  while (left < right) {
    size_t mid = left + (right - left) / 2;
    SimilarityRecord *record = &records->items[mid];
    int cmp = strcmp(record->query_id, query_id);
    if (cmp == 0) {
      cmp = strcmp(record->subject_id, subject_id);
    }
    if (cmp < 0) {
      left = mid + 1;
    } else {
      right = mid;
    }
  }
  if (left < records->len &&
      strcmp(records->items[left].query_id, query_id) == 0 &&
      strcmp(records->items[left].subject_id, subject_id) == 0) {
    return &records->items[left];
  }
  return NULL;
}

static bool has_string(StringVec *vec, const char *value) {
  for (size_t i = 0; i < vec->len; i++) {
    if (strcmp(vec->items[i], value) == 0) {
      return true;
    }
  }
  return false;
}

static char *dup_orient(const char *value) {
  return xstrdup(value);
}

static void free_edges(EdgeVec *edges) {
  for (size_t i = 0; i < edges->len; i++) {
    free(edges->items[i].seq_a);
    free(edges->items[i].seq_b);
    free(edges->items[i].taxon_a);
    free(edges->items[i].taxon_b);
  }
  free(edges->items);
}

static void build_orthologs(RecordVec *records, double percent_match_cutoff, int evalue_exp_cutoff, EdgeVec *out_edges) {
  for (size_t i = 0; i < records->len; i++) {
    SimilarityRecord *left = &records->items[i];
    if (strcmp(left->query_taxon, left->subject_taxon) == 0) continue;
    if (!passes_thresholds(left, percent_match_cutoff, evalue_exp_cutoff)) continue;

    SimilarityRecord *best = NULL;
    for (size_t j = 0; j < records->len; j++) {
      SimilarityRecord *candidate = &records->items[j];
      if (strcmp(candidate->query_id, left->query_id) != 0) continue;
      if (strcmp(candidate->subject_taxon, left->subject_taxon) != 0) continue;
      if (strcmp(candidate->query_taxon, candidate->subject_taxon) == 0) continue;
      if (!best || compare_evalue(candidate, best) < 0) {
        best = candidate;
      }
    }
    if (!best) continue;
    if (!(left->evalue_mant < 0.01 || compare_evalue(left, best) == 0)) continue;
    if (strcmp(left->query_id, left->subject_id) >= 0) continue;

    SimilarityRecord *right = find_pair(records, left->subject_id, left->query_id);
    if (!right) continue;
    if (strcmp(right->query_taxon, right->subject_taxon) == 0) continue;
    if (!passes_thresholds(right, percent_match_cutoff, evalue_exp_cutoff)) continue;

    SimilarityRecord *reverse_best = NULL;
    for (size_t j = 0; j < records->len; j++) {
      SimilarityRecord *candidate = &records->items[j];
      if (strcmp(candidate->query_id, right->query_id) != 0) continue;
      if (strcmp(candidate->subject_taxon, right->subject_taxon) != 0) continue;
      if (strcmp(candidate->query_taxon, candidate->subject_taxon) == 0) continue;
      if (!reverse_best || compare_evalue(candidate, reverse_best) < 0) {
        reverse_best = candidate;
      }
    }
    if (!reverse_best) continue;
    if (!(right->evalue_mant < 0.01 || compare_evalue(right, reverse_best) == 0)) continue;

    EdgeRecord edge;
    edge.seq_a = dup_orient(left->query_id);
    edge.seq_b = dup_orient(left->subject_id);
    edge.taxon_a = dup_orient(left->query_taxon);
    edge.taxon_b = dup_orient(left->subject_taxon);
    edge.unnormalized_score = score_from_records(left, right, 0.01);
    edge.normalized_score = 0.0;
    push_edge(out_edges, edge);
  }

  for (size_t i = 0; i < out_edges->len; i++) {
    DoubleVec scores = {0};
    EdgeRecord *edge = &out_edges->items[i];
    for (size_t j = 0; j < out_edges->len; j++) {
      EdgeRecord *candidate = &out_edges->items[j];
      bool same_pair =
          (strcmp(edge->taxon_a, candidate->taxon_a) == 0 && strcmp(edge->taxon_b, candidate->taxon_b) == 0) ||
          (strcmp(edge->taxon_a, candidate->taxon_b) == 0 && strcmp(edge->taxon_b, candidate->taxon_a) == 0);
      if (same_pair) {
        push_double(&scores, candidate->unnormalized_score);
      }
    }
    double sum = 0.0;
    for (size_t k = 0; k < scores.len; k++) sum += scores.items[k];
    edge->normalized_score = edge->unnormalized_score / (sum / (double)scores.len);
    free(scores.items);
  }

  qsort(out_edges->items, out_edges->len, sizeof(EdgeRecord), cmp_edge_output);
}

static void build_inparalogs(RecordVec *records, EdgeVec *orthologs, double percent_match_cutoff, int evalue_exp_cutoff, EdgeVec *out_edges) {
  StringVec ortholog_ids = {0};
  for (size_t i = 0; i < orthologs->len; i++) {
    (void)string_index(&ortholog_ids, orthologs->items[i].seq_a);
    (void)string_index(&ortholog_ids, orthologs->items[i].seq_b);
  }

  for (size_t i = 0; i < records->len; i++) {
    SimilarityRecord *left = &records->items[i];
    if (strcmp(left->query_taxon, left->subject_taxon) != 0) continue;
    if (strcmp(left->query_id, left->subject_id) == 0) continue;
    if (!passes_thresholds(left, percent_match_cutoff, evalue_exp_cutoff)) continue;

    SimilarityRecord *best_inter = NULL;
    for (size_t j = 0; j < records->len; j++) {
      SimilarityRecord *candidate = &records->items[j];
      if (strcmp(candidate->query_id, left->query_id) != 0) continue;
      if (strcmp(candidate->query_taxon, candidate->subject_taxon) == 0) continue;
      if (!best_inter || compare_evalue(candidate, best_inter) < 0) {
        best_inter = candidate;
      }
    }
    if (best_inter) {
      bool keep = left->evalue_mant < 0.001 ||
                  left->evalue_exp < best_inter->evalue_exp ||
                  (left->evalue_exp == best_inter->evalue_exp && left->evalue_mant <= best_inter->evalue_mant);
      if (!keep) continue;
    }
    if (strcmp(left->query_id, left->subject_id) >= 0) continue;

    SimilarityRecord *right = find_pair(records, left->subject_id, left->query_id);
    if (!right) continue;
    if (strcmp(right->query_taxon, right->subject_taxon) != 0) continue;
    if (!passes_thresholds(right, percent_match_cutoff, evalue_exp_cutoff)) continue;

    SimilarityRecord *right_best_inter = NULL;
    for (size_t j = 0; j < records->len; j++) {
      SimilarityRecord *candidate = &records->items[j];
      if (strcmp(candidate->query_id, right->query_id) != 0) continue;
      if (strcmp(candidate->query_taxon, candidate->subject_taxon) == 0) continue;
      if (!right_best_inter || compare_evalue(candidate, right_best_inter) < 0) {
        right_best_inter = candidate;
      }
    }
    if (right_best_inter) {
      bool keep = right->evalue_mant < 0.001 ||
                  right->evalue_exp < right_best_inter->evalue_exp ||
                  (right->evalue_exp == right_best_inter->evalue_exp && right->evalue_mant <= right_best_inter->evalue_mant);
      if (!keep) continue;
    }

    EdgeRecord edge;
    edge.seq_a = dup_orient(left->query_id);
    edge.seq_b = dup_orient(left->subject_id);
    edge.taxon_a = dup_orient(left->query_taxon);
    edge.taxon_b = dup_orient(left->query_taxon);
    edge.unnormalized_score = score_from_records(left, right, 0.01);
    edge.normalized_score = 0.0;
    push_edge(out_edges, edge);
  }

  for (size_t i = 0; i < out_edges->len; i++) {
    EdgeRecord *edge = &out_edges->items[i];
    DoubleVec all_scores = {0};
    DoubleVec orth_scores = {0};
    for (size_t j = 0; j < out_edges->len; j++) {
      EdgeRecord *candidate = &out_edges->items[j];
      if (strcmp(edge->taxon_a, candidate->taxon_a) != 0) continue;
      push_double(&all_scores, candidate->unnormalized_score);
      if (has_string(&ortholog_ids, candidate->seq_a) || has_string(&ortholog_ids, candidate->seq_b)) {
        push_double(&orth_scores, candidate->unnormalized_score);
      }
    }
    DoubleVec *chosen = orth_scores.len ? &orth_scores : &all_scores;
    double sum = 0.0;
    for (size_t k = 0; k < chosen->len; k++) sum += chosen->items[k];
    edge->normalized_score = edge->unnormalized_score / (sum / (double)chosen->len);
    free(all_scores.items);
    free(orth_scores.items);
  }

  qsort(out_edges->items, out_edges->len, sizeof(EdgeRecord), cmp_edge_output);

  for (size_t i = 0; i < ortholog_ids.len; i++) free(ortholog_ids.items[i]);
  free(ortholog_ids.items);
}

static bool edge_exists(EdgeVec *edges, const char *seq_a, const char *seq_b) {
  for (size_t i = 0; i < edges->len; i++) {
    if (strcmp(edges->items[i].seq_a, seq_a) == 0 && strcmp(edges->items[i].seq_b, seq_b) == 0) return true;
  }
  return false;
}

static void maybe_add_coortholog_candidate(
    RecordVec *records,
    EdgeVec *orthologs,
    EdgeVec *out_edges,
    const char *seq_a,
    const char *seq_b,
    double percent_match_cutoff,
    int evalue_exp_cutoff
) {
  const char *first = strcmp(seq_a, seq_b) < 0 ? seq_a : seq_b;
  const char *second = strcmp(seq_a, seq_b) < 0 ? seq_b : seq_a;
  if (strcmp(first, second) == 0) return;
  if (edge_exists(orthologs, first, second)) return;
  if (edge_exists(out_edges, first, second)) return;

  SimilarityRecord *left = find_pair(records, first, second);
  SimilarityRecord *right = find_pair(records, second, first);
  if (!left || !right) return;
  if (!passes_thresholds(left, percent_match_cutoff, evalue_exp_cutoff)) return;
  if (!passes_thresholds(right, percent_match_cutoff, evalue_exp_cutoff)) return;

  EdgeRecord edge;
  edge.seq_a = dup_orient(first);
  edge.seq_b = dup_orient(second);
  edge.taxon_a = dup_orient(left->query_taxon);
  edge.taxon_b = dup_orient(left->subject_taxon);
  edge.unnormalized_score = score_from_records(left, right, 0.00001);
  edge.normalized_score = 0.0;
  push_edge(out_edges, edge);
}

static void build_coorthologs(RecordVec *records, EdgeVec *orthologs, EdgeVec *inparalogs, double percent_match_cutoff, int evalue_exp_cutoff, EdgeVec *out_edges) {
  for (size_t i = 0; i < inparalogs->len; i++) {
    const char *ip_left[2] = {inparalogs->items[i].seq_a, inparalogs->items[i].seq_b};
    const char *ip_right[2] = {inparalogs->items[i].seq_b, inparalogs->items[i].seq_a};
    for (size_t j = 0; j < orthologs->len; j++) {
      const char *ortho_left[2] = {orthologs->items[j].seq_a, orthologs->items[j].seq_b};
      const char *ortho_right[2] = {orthologs->items[j].seq_b, orthologs->items[j].seq_a};
      for (int idir = 0; idir < 2; idir++) {
        for (int odir = 0; odir < 2; odir++) {
          if (strcmp(ip_right[idir], ortho_left[odir]) == 0) {
            maybe_add_coortholog_candidate(
                records,
                orthologs,
                out_edges,
                ip_left[idir],
                ortho_right[odir],
                percent_match_cutoff,
                evalue_exp_cutoff);
          }
        }
      }
    }
  }

  for (size_t i = 0; i < inparalogs->len; i++) {
    const char *ip1_left[2] = {inparalogs->items[i].seq_a, inparalogs->items[i].seq_b};
    const char *ip1_right[2] = {inparalogs->items[i].seq_b, inparalogs->items[i].seq_a};
    for (size_t j = 0; j < orthologs->len; j++) {
      const char *ortho_left[2] = {orthologs->items[j].seq_a, orthologs->items[j].seq_b};
      const char *ortho_right[2] = {orthologs->items[j].seq_b, orthologs->items[j].seq_a};
      for (size_t k = 0; k < inparalogs->len; k++) {
        const char *ip2_left[2] = {inparalogs->items[k].seq_a, inparalogs->items[k].seq_b};
        const char *ip2_right[2] = {inparalogs->items[k].seq_b, inparalogs->items[k].seq_a};
        for (int idir1 = 0; idir1 < 2; idir1++) {
          for (int odir = 0; odir < 2; odir++) {
            if (strcmp(ip1_right[idir1], ortho_left[odir]) != 0) continue;
            for (int idir2 = 0; idir2 < 2; idir2++) {
              if (strcmp(ortho_right[odir], ip2_left[idir2]) == 0) {
                maybe_add_coortholog_candidate(
                    records,
                    orthologs,
                    out_edges,
                    ip1_left[idir1],
                    ip2_right[idir2],
                    percent_match_cutoff,
                    evalue_exp_cutoff);
              }
            }
          }
        }
      }
    }
  }

  for (size_t i = 0; i < out_edges->len; i++) {
    DoubleVec scores = {0};
    EdgeRecord *edge = &out_edges->items[i];
    for (size_t j = 0; j < out_edges->len; j++) {
      EdgeRecord *candidate = &out_edges->items[j];
      bool same_pair =
          (strcmp(edge->taxon_a, candidate->taxon_a) == 0 && strcmp(edge->taxon_b, candidate->taxon_b) == 0) ||
          (strcmp(edge->taxon_a, candidate->taxon_b) == 0 && strcmp(edge->taxon_b, candidate->taxon_a) == 0);
      if (same_pair) push_double(&scores, candidate->unnormalized_score);
    }
    double sum = 0.0;
    for (size_t k = 0; k < scores.len; k++) sum += scores.items[k];
    edge->normalized_score = edge->unnormalized_score / (sum / (double)scores.len);
    free(scores.items);
  }

  qsort(out_edges->items, out_edges->len, sizeof(EdgeRecord), cmp_edge_output);
}

static void write_edges(const char *path, EdgeVec *edges) {
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open output file");
  for (size_t i = 0; i < edges->len; i++) {
    fprintf(handle, "%s\t%s\t%.3f\n", edges->items[i].seq_a, edges->items[i].seq_b, rounded_score(edges->items[i].normalized_score));
  }
  fclose(handle);
}

static void write_mcl_input(const char *path, EdgeVec *orthologs, EdgeVec *inparalogs, EdgeVec *coorthologs) {
  EdgeVec merged = {0};
  for (size_t i = 0; i < orthologs->len; i++) push_edge(&merged, orthologs->items[i]);
  for (size_t i = 0; i < inparalogs->len; i++) push_edge(&merged, inparalogs->items[i]);
  for (size_t i = 0; i < coorthologs->len; i++) push_edge(&merged, coorthologs->items[i]);
  qsort(merged.items, merged.len, sizeof(EdgeRecord), cmp_edge_union);
  FILE *handle = fopen(path, "w");
  if (!handle) die("Could not open mclInput");
  for (size_t i = 0; i < merged.len; i++) {
    fprintf(handle, "%s\t%s\t%.3f\n", merged.items[i].seq_a, merged.items[i].seq_b, rounded_score(merged.items[i].normalized_score));
  }
  fclose(handle);
  free(merged.items);
}

static void make_dir(const char *path) {
  char command[4096];
  snprintf(command, sizeof(command), "mkdir -p \"%s\"", path);
  int rc = system(command);
  if (rc != 0) {
    die("Could not create output directory");
  }
}

static RecordVec read_records(const char *path) {
  RecordVec records = {0};
  FILE *handle = fopen(path, "r");
  if (!handle) die("Could not open similarSequences file");
  char *line = NULL;
  size_t cap = 0;
  while (getline(&line, &cap, handle) != -1) {
    if (line[0] == '\n' || line[0] == '\0') continue;
    char *fields[8];
    int idx = 0;
    char *saveptr = NULL;
    char *token = strtok_r(line, "\t\n", &saveptr);
    while (token && idx < 8) {
      fields[idx++] = token;
      token = strtok_r(NULL, "\t\n", &saveptr);
    }
    if (idx != 8) die("similarSequences line does not have 8 columns");
    SimilarityRecord record;
    record.query_id = xstrdup(fields[0]);
    record.subject_id = xstrdup(fields[1]);
    record.query_taxon = xstrdup(fields[2]);
    record.subject_taxon = xstrdup(fields[3]);
    record.evalue_mant = atof(fields[4]);
    record.evalue_exp = atoi(fields[5]);
    record.percent_identity = atof(fields[6]);
    record.percent_match = atof(fields[7]);
    push_record(&records, record);
  }
  free(line);
  fclose(handle);
  qsort(records.items, records.len, sizeof(SimilarityRecord), cmp_record_pair);
  return records;
}

static void free_records(RecordVec *records) {
  for (size_t i = 0; i < records->len; i++) {
    free(records->items[i].query_id);
    free(records->items[i].subject_id);
    free(records->items[i].query_taxon);
    free(records->items[i].subject_taxon);
  }
  free(records->items);
}

static void usage(void) {
  fprintf(stderr, "usage: pairs_engine similar_sequences out_dir percent_match_cutoff evalue_exp_cutoff\n");
}

int main(int argc, char **argv) {
  if (argc != 5) {
    usage();
    return 1;
  }

  const char *similar_sequences = argv[1];
  const char *out_dir = argv[2];
  double percent_match_cutoff = atof(argv[3]);
  int evalue_exp_cutoff = atoi(argv[4]);

  RecordVec records = read_records(similar_sequences);
  EdgeVec orthologs = {0};
  EdgeVec inparalogs = {0};
  EdgeVec coorthologs = {0};

  build_orthologs(&records, percent_match_cutoff, evalue_exp_cutoff, &orthologs);
  build_inparalogs(&records, &orthologs, percent_match_cutoff, evalue_exp_cutoff, &inparalogs);
  build_coorthologs(&records, &orthologs, &inparalogs, percent_match_cutoff, evalue_exp_cutoff, &coorthologs);

  char pairs_dir[4096];
  snprintf(pairs_dir, sizeof(pairs_dir), "%s/pairs", out_dir);
  make_dir(out_dir);
  make_dir(pairs_dir);

  char orthologs_path[4096];
  char inparalogs_path[4096];
  char coorthologs_path[4096];
  char mcl_input_path[4096];
  snprintf(orthologs_path, sizeof(orthologs_path), "%s/orthologs.txt", pairs_dir);
  snprintf(inparalogs_path, sizeof(inparalogs_path), "%s/inparalogs.txt", pairs_dir);
  snprintf(coorthologs_path, sizeof(coorthologs_path), "%s/coorthologs.txt", pairs_dir);
  snprintf(mcl_input_path, sizeof(mcl_input_path), "%s/mclInput", out_dir);

  write_edges(orthologs_path, &orthologs);
  write_edges(inparalogs_path, &inparalogs);
  write_edges(coorthologs_path, &coorthologs);
  write_mcl_input(mcl_input_path, &orthologs, &inparalogs, &coorthologs);

  free_edges(&orthologs);
  free_edges(&inparalogs);
  free_edges(&coorthologs);
  free_records(&records);
  return 0;
}
