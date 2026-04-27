#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>

// macros
#define SNP_NUM_COLS 6
#define SNP_VAR_COL 0
#define SNP_CHR_COL 1
#define SNP_CM_COL 2
#define SNP_POS_COL 3
#define SNP_REF_COL 4
#define SNP_ALT_COL 5

#define IND_NUM_COLS 3
#define IND_ID_COL 0
#define IND_SEX_COL 1
#define IND_POP_COL 2

#define BITS_IN_BYTE 8
#define RECORD_ELEM_SIZE_BITS 2
#define RECORD_ELEMS_PER_BYTE 4
#define RECORD_ELEMS_MASK_BASE 3

/* BASIC TYPES */
typedef struct {
	size_t length;
	char** var_id;
	char** chr;
	double* cm;
	uint64_t* pos;
	char** ref;
	char** alt;
	uint32_t hash;
} snp_data;

typedef struct {
	size_t length;
	char** ind_id;
	char** sex;
	char** population;
	uint32_t hash;
} ind_data;

typedef struct {
	bool is_hdr_read;
	size_t idx;
	size_t record_size;
	FILE* fp;
} pam_file;  // PACKEDANCESTRYMAP reader

typedef struct {
	size_t idx;
	size_t record_size;
	FILE* fp;
} egn_file;  // EIGENSTRAT reader

typedef struct {
	size_t n_ind;
	size_t n_snp;
	uint32_t ind_hash;
	uint32_t snp_hash;
} hdr_data;

/* BASIC FUNCTIONS */
uint32_t hash_str(char* str) {
	uint32_t hash_out = 0;
	size_t str_len = strlen(str);
	for(int i = 0; i < str_len; i++) {
		hash_out *= 23;
		hash_out += str[i];
	}
	return(hash_out);
}

FILE* safe_read(char* filename, char* mode) {
	FILE *fp = fopen(filename, mode);
	if (fp == NULL) {
		fprintf(stderr, "ERROR: cannot open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
	return fp;
}

size_t get_filesize(char* filename) {
	FILE *fp = safe_read(filename, "rb");
	fseek(fp, 0, SEEK_END);
	size_t file_size = ftell(fp);
	fclose(fp);
	return file_size;
}

uint64_t get_number_of_lines(char* filename) {
	uint64_t num_lines = 0;
	FILE *fp = safe_read(filename, "r");
	while(!feof(fp)) {
		num_lines += (fgetc(fp) == '\n');
	}
	fclose(fp);
	return num_lines;
}

char** get_column_elems(char* col_str, size_t ncol) {
	char* cur_col = strtok(col_str, " \n\t\r");
	size_t col_num = 0;
	char** col_elems = (char**) malloc(ncol * sizeof(char*));
	while(cur_col != NULL) {
		if(col_num == ncol) {
			if(strcmp(cur_col, "") == 0) {
				return col_elems;
			} else {
				fprintf(stderr, "ERROR: Too many columns in '.snp' or '.ind' file.\n");
				exit(EXIT_FAILURE);
			}
		} else {
			col_elems[col_num] = cur_col;
			col_num += 1;
			cur_col = strtok(NULL, " \n\t\r");
		}
	}
	if (col_num == ncol) {
		return col_elems;
	} else {
		fprintf(stderr, "ERROR: Too few columns in '.snp' or '.ind' file.\n");
		exit(EXIT_FAILURE);
	}
}

snp_data read_snp_file(char* filename) {
	FILE *fp = safe_read(filename, "r");
	uint64_t num_snps = get_number_of_lines(filename);
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	
	snp_data snp_info = {
		.length = num_snps,
		.var_id = (char**) malloc(num_snps * sizeof(char*)),
		.chr = (char**) malloc(num_snps * sizeof(char*)),
		.cm = (double*) malloc(num_snps * sizeof(double)),
		.pos = (uint64_t*) malloc(num_snps * sizeof(uint64_t)),
		.ref = (char**) malloc(num_snps * sizeof(char*)),
		.alt = (char**) malloc(num_snps * sizeof(char*)),
		.hash = 0
	};

	size_t idx = 0;
	while((nread = getline(&line, &size, fp)) != -1) {
		char** elems = get_column_elems(line, SNP_NUM_COLS);
		snp_info.var_id[idx] = strdup(elems[SNP_VAR_COL]);
		snp_info.chr[idx] = strdup(elems[SNP_CHR_COL]);
		sscanf(elems[SNP_CM_COL], "%lf", &snp_info.cm[idx]);
		snp_info.pos[idx] = atoi(elems[SNP_POS_COL]);
		snp_info.ref[idx] = strdup(elems[SNP_REF_COL]);
		snp_info.alt[idx] = strdup(elems[SNP_ALT_COL]);
		// calculate hash
		snp_info.hash *= 17;
		snp_info.hash = snp_info.hash ^ hash_str(snp_info.var_id[idx]);
		idx++;
	}
	free(line);
	fclose(fp);
	return snp_info;
}

ind_data read_ind_file(char* filename) {
	FILE *fp = safe_read(filename, "r");
	uint64_t num_inds = get_number_of_lines(filename);
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	
	ind_data ind_info = {
		.length = num_inds,
		.ind_id = (char**) malloc(num_inds * sizeof(char*)),
		.sex = (char**) malloc(num_inds * sizeof(char*)),
		.population = (char**) malloc(num_inds * sizeof(char*)),
		.hash = 0
	};

	size_t idx = 0;
	while((nread = getline(&line, &size, fp)) != -1) {
		char** elems = get_column_elems(line, IND_NUM_COLS);
		ind_info.ind_id[idx] = strdup(elems[IND_ID_COL]);
		ind_info.sex[idx] = strdup(elems[IND_SEX_COL]);
		ind_info.population[idx] = strdup(elems[IND_POP_COL]);
		// calculate hash
		ind_info.hash *= 17;
		ind_info.hash = ind_info.hash ^ hash_str(ind_info.ind_id[idx]);
		idx++;
	}
	free(line);
	fclose(fp);
	return ind_info;
}

pam_file open_pam(char* filename, size_t n_snp) {
	pam_file pf;
	size_t file_size = get_filesize(filename);
	pf.is_hdr_read = false;
	pf.idx = 0;
	pf.record_size = file_size / (n_snp + 1);
	if ((pf.record_size * (n_snp+1)) != file_size) {
		fprintf(stderr, "Invalid PACKEDANCESTRYMAP file. File size must be a multiple of the number of SNPs.\n");
		exit(EXIT_FAILURE);
	}
	pf.fp = safe_read(filename, "rb");
	return pf;
}

hdr_data read_pam_header(pam_file* pf) {
	if(pf->is_hdr_read) {
		fprintf(stderr, "Header already read!");
		exit(EXIT_FAILURE);
	}
	hdr_data hdr_info;
	char* hdr = (char*) malloc(sizeof(char) * pf->record_size);
	fread(hdr, 1, pf->record_size, pf->fp);
	sscanf(hdr, "GENO   %u %u %x %x", &hdr_info.n_ind, &hdr_info.n_snp, &hdr_info.ind_hash, &hdr_info.snp_hash);
	free(hdr);
	pf->is_hdr_read = true;
	return hdr_info;
}

uint8_t* read_pam_record(pam_file* pf, size_t n_ind) {
	if(!pf->is_hdr_read) {
		fprintf(stderr, "Header must be read before reading records!\n");
		exit(EXIT_FAILURE);
	}
	uint8_t* record = (uint8_t*)malloc(n_ind * sizeof(uint8_t));
	size_t num_leftover_bytes = pf->record_size - (int)ceil((float)(n_ind*RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE);
	uint8_t record_byte;
	for(int i = 0; i < n_ind; i++) {
		uint8_t elem_pos = i%RECORD_ELEMS_PER_BYTE;
		uint8_t shift_by = (BITS_IN_BYTE-RECORD_ELEM_SIZE_BITS) - (RECORD_ELEM_SIZE_BITS*elem_pos);
		if(elem_pos == 0) {
			record_byte = getc(pf->fp);
		}
		record[i] = (record_byte & (RECORD_ELEMS_MASK_BASE << shift_by)) >> shift_by;
	}
	fseek(pf->fp, num_leftover_bytes, SEEK_CUR);  // skip useless bytes
	return record;
}

bool check_ind_hash(ind_data* ind_info, hdr_data* hdr_info) {
	return ind_info->hash == hdr_info->ind_hash;
}

bool check_snp_hash(snp_data* snp_info, hdr_data* hdr_info) {
	return snp_info->hash == hdr_info->ind_hash;
}
