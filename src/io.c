#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

// macros
#define SNP_NUM_COLS 6
#define SNP_VAR_COL 0
#define SNP_CHR_COL 1
#define SNP_CM_COL 2
#define SNP_POS_COL 3
#define SNP_REF_COL 4
#define SNP_ALT_COL 5

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
	size_t start;
	size_t length;
} col_pos;

uint32_t hash_str(char* str) {
	uint32_t hash_out = 0;
	for(int i = 0; i < strlen(str); i++) {
		hash_out *= 23;
		hash_out += str[i];
	}
	return(hash_out);
}

size_t get_filesize(char* filename) {
	FILE *fp = fopen(filename, "rb");
	fseek(fp, 0, SEEK_END);
	size_t file_size = ftell(fp);
	fclose(fp);
	return file_size;
}

uint64_t get_number_of_lines(char* filename) {
	uint64_t num_lines = 0;
	FILE *fp = fopen(filename, "r");
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
				printf("ERROR: Too many columns in '.snp' or '.ind' file.\n");
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
		printf("ERROR: Too few columns in '.snp' or '.ind' file.\n"); exit(EXIT_FAILURE);
	}
}

snp_data read_snp_file(char* filename) {
	FILE *fp = fopen(filename, "r");
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
