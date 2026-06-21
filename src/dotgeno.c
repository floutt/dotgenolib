#define  _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include <sys/queue.h>
#include "khash.h"
#include "dotgeno.h"

// macros
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

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

/* hash table stuff */
/**
 *  @brief Hash function for ind_idx for use with khash.h, based on the djb2 hash function
 *
 *  @param[in] ind_idx_struct ind_idx structure to be hashed
 *
 *  @return hash value of type khint_t
 */
static inline khint_t hash_ind_idx(ind_idx ind_idx_struct) {
	// String linker for .ind hashing
	int IND_LINK_LEN = 20;
	char IND_LINK[21] = "gzvrEy55bcEN0gqRqvL6";

	khint_t hash = 5381;
	
	// first string
	int len = strlen(ind_idx_struct.ind_id);
	for(int i = 0; i < len; i++) {
		hash = ((hash << 5) + hash) + ind_idx_struct.ind_id[i];
	}
	
	// linker string
	for(int i = 0; i < IND_LINK_LEN; i++) {
		hash = ((hash << 5) + hash) + IND_LINK[i];
	}

	// second string
	len = strlen(ind_idx_struct.ind_pop);
	for(int i = 0; i < len; i++) {
		hash = ((hash << 5) + hash) + ind_idx_struct.ind_pop[i];
	}
	return hash;
}

/*! @brief comparison function for ind_idx hash. For use with khash.h */
#define ind_idx_equal(a, b) ((strcmp((a).ind_id, (b).ind_id) == 0) && strcmp((a).ind_pop, (b).ind_pop) == 0)

KHASH_MAP_INIT_STR(ID_MAP_STR, size_t)
KHASH_INIT(ID_MAP_IND, ind_idx, size_t, true, hash_ind_idx, ind_idx_equal)

struct id_map_str {
    khash_t(ID_MAP_STR)* map;
};

struct id_map_ind {
    khash_t(ID_MAP_IND)* map;
};

/* BASIC IO FUNCTIONS */
FILE* safe_read(char* filename, char* mode) {
	FILE* fp = fopen(filename, mode);
	if (fp == NULL) {
		fprintf(stderr, "ERROR: cannot open file %s\n", filename);
		exit(EXIT_FAILURE);
	}
	return fp;
}

size_t get_filesize(char* filename) {
	FILE* fp = safe_read(filename, "rb");
	fseek(fp, 0, SEEK_END);
	size_t file_size = ftell(fp);
	fclose(fp);
	return file_size;
}

uint64_t get_number_of_lines(char* filename) {
	uint64_t num_lines = 0;
	FILE* fp = safe_read(filename, "r");
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

/* Specific EIGENSTRAT and PACKEDANCESTRYMAP functions */
/**
 *  @brief Hash function for single string for use in evaluating integrity of .ind and .snp files
 *
 *  This function is part of the hashing process used for these files. This one operates on
 *  individual IDs or variant IDs. The hash is built as the .snp and .ind files are being read.
 *  See read_snp_file and read_ind_file for more details on their implementation.
 *
 *  @param[in] str string to be hashed
 *
 *  @return hash value of type uint32_t
 */
static inline uint32_t hash_str(char* str) {
	uint32_t hash_out = 0;
	size_t str_len = strlen(str);
	for(size_t i = 0; i < str_len; i++) {
		hash_out *= 23;
		hash_out += str[i];
	}
	return(hash_out);
}

snp_data read_snp_file(char* filename) {
	FILE* fp = safe_read(filename, "r");
	uint64_t num_snps = get_number_of_lines(filename);
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	
	snp_data snp_info = {
		.length = num_snps,
		.var_id = (char**) malloc(num_snps * sizeof(char*)),
		.chr = (char**) malloc(num_snps * sizeof(char*)),
		.pos_morgans = (double*) malloc(num_snps * sizeof(double)),
		.pos = (uint64_t*) malloc(num_snps * sizeof(uint64_t)),
		.ref = (char**) malloc(num_snps * sizeof(char*)),
		.alt = (char**) malloc(num_snps * sizeof(char*)),
		.rev_idx = (id_map_str*) malloc(sizeof(id_map_str)),
		.hash = 0
	};
	snp_info.rev_idx->map = kh_init(ID_MAP_STR);
	size_t idx = 0;
	while((nread = getline(&line, &size, fp)) != -1) {
		char** elems = get_column_elems(line, SNP_NUM_COLS);
		snp_info.var_id[idx] = strdup(elems[SNP_VAR_COL]);
		snp_info.chr[idx] = strdup(elems[SNP_CHR_COL]);
		sscanf(elems[SNP_CM_COL], "%lf", &snp_info.pos_morgans[idx]);
		snp_info.pos[idx] = atoi(elems[SNP_POS_COL]);
		snp_info.ref[idx] = strdup(elems[SNP_REF_COL]);
		snp_info.alt[idx] = strdup(elems[SNP_ALT_COL]);
		// calculate hash
		snp_info.hash *= 17;
		snp_info.hash = snp_info.hash ^ hash_str(snp_info.var_id[idx]);
		free(elems);
		int ret;
		khiter_t k = kh_put(ID_MAP_STR, snp_info.rev_idx->map, snp_info.var_id[idx], &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert variant name into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Multiple instances of variant %s found in the .snp file.\n", snp_info.var_id[idx]);
			exit(EXIT_FAILURE);
		}
		kh_value(snp_info.rev_idx->map, k) = idx;
		idx++;
	}
	free(line);
	fclose(fp);
	return snp_info;
}

ind_data read_ind_file(char* filename) {
	FILE* fp = safe_read(filename, "r");
	uint64_t num_inds = get_number_of_lines(filename);
	char* line = NULL;
	size_t size = 0;
	ssize_t nread;
	
	ind_data ind_info = {
		.length = num_inds,
		.ind_id = (char**) malloc(num_inds * sizeof(char*)),
		.sex = (char**) malloc(num_inds * sizeof(char*)),
		.population = (char**) malloc(num_inds * sizeof(char*)),
		.rev_idx = (id_map_ind*) malloc(sizeof(id_map_ind)),
		.hash = 0
	};
	ind_info.rev_idx->map = kh_init(ID_MAP_IND);
	size_t idx = 0;
	while((nread = getline(&line, &size, fp)) != -1) {
		char** elems = get_column_elems(line, IND_NUM_COLS);
		ind_info.ind_id[idx] = strdup(elems[IND_ID_COL]);
		ind_info.sex[idx] = strdup(elems[IND_SEX_COL]);
		ind_info.population[idx] = strdup(elems[IND_POP_COL]);
		free(elems);
		// calculate hash
		ind_info.hash *= 17;
		ind_info.hash = ind_info.hash ^ hash_str(ind_info.ind_id[idx]);
		int ret;
		ind_idx ind_idx_struct = (ind_idx){ind_info.ind_id[idx], ind_info.population[idx]};
		khiter_t k = kh_put(ID_MAP_IND, ind_info.rev_idx->map, ind_idx_struct, &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert individual into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Multiple instances of individual (%s, %s) found in the .ind file.\n", ind_info.ind_id[idx], ind_info.population[idx]);
			exit(EXIT_FAILURE);
		}
		kh_value(ind_info.rev_idx->map, k) = idx;
		idx++;
	}
	free(line);
	fclose(fp);
	return ind_info;
}

short write_snp_data(snp_data* snp_info, char* filename) {
	FILE* fp = safe_read(filename, "w+");
	for(size_t i = 0; i < snp_info->length; i++) {
		fprintf(fp, "%s\t%s\t%lf\t%zu\t%s\t%s\n",
				snp_info->var_id[i], snp_info->chr[i], snp_info->pos_morgans[i], snp_info->pos[i], snp_info->ref[i], snp_info->alt[i]);
	}
	return fclose(fp);
}

short write_ind_data(ind_data* ind_info, char* filename) {
	FILE* fp = safe_read(filename, "w+");
	for(size_t i = 0; i < ind_info->length; i++) {
		fprintf(fp, "%s\t%s\t%s\n",
				ind_info->ind_id[i], ind_info->sex[i], ind_info->population[i]);
	}
	return fclose(fp);
}

short get_snp_idx(snp_data* snp_info, char* var_name, size_t* idx) {
	khint_t k = kh_get(ID_MAP_STR, snp_info->rev_idx->map, var_name);
	if(k == kh_end(snp_info->rev_idx->map)) {
		return -1;
	} else {
		*idx = kh_value(snp_info->rev_idx->map, k);
		return 0;
	}
}

short get_ind_idx(ind_data* ind_info, char* ind_id, char* ind_pop, size_t* idx) {
	ind_idx iidx;
	iidx.ind_id = ind_id;
	iidx.ind_pop = ind_pop; 
	khint_t k = kh_get(ID_MAP_IND, ind_info->rev_idx->map, iidx);
	if(k == kh_end(ind_info->rev_idx->map)) {
		return -1;
	} else {
		*idx = kh_value(ind_info->rev_idx->map, k);
		return 0;
	}
}

void get_multiple_snp_idx(snp_data* snp_info, char** var_names, size_t length, struct idx_head* head_idx, struct str_list_head* head_str) {
	for(size_t i = 0; i < length; i++) {
		size_t idx;
		short ret = get_snp_idx(snp_info, var_names[i], &idx);
		if(ret == -1) {
			if(head_str) {
				struct str_node* stn = (struct str_node*)malloc(sizeof(struct str_node));
				stn->str = strdup(var_names[i]);
				STAILQ_INSERT_TAIL(head_str, stn, nodes);
			}
			continue;
		}
		struct idx_node* idn = (struct idx_node*)malloc(sizeof(struct idx_node));
		idn->idx = idx;
		STAILQ_INSERT_TAIL(head_idx, idn, nodes);
	}
}

void get_multiple_ind_idx(ind_data* ind_info, char** ind_ids, char** ind_pops, size_t length, struct idx_head* head_idx, struct ind_idx_head* head_iidx) {
	for(size_t i = 0; i < length; i++) {
		size_t idx;
		short ret = get_ind_idx(ind_info, ind_ids[i], ind_pops[i], &idx);
		if(ret == -1) {
			if(head_iidx) {
				struct ind_idx_node* idn = (struct ind_idx_node*)malloc(sizeof(struct ind_idx_node));
				idn->iidx = (ind_idx*)malloc(sizeof(ind_idx));
				idn->iidx->ind_id = strdup(ind_ids[i]);
				idn->iidx->ind_pop = strdup(ind_pops[i]);
				STAILQ_INSERT_TAIL(head_iidx, idn, nodes);
			}
			continue;
		}
		struct idx_node* idn = (struct idx_node*)malloc(sizeof(struct idx_node));
		idn->idx = idx;
		STAILQ_INSERT_TAIL(head_idx, idn, nodes);
	}
}

short filter_snp_data(snp_data* snp_in, snp_data* snp_out, struct idx_head* head) {
	// get length
	size_t length = 0;
	struct idx_node* tmp_node;
	STAILQ_FOREACH(tmp_node, head, nodes) {
		length++;
	}
	if(length == 0) {
		return -1;
	}
	// make struct
	snp_out->length = length,
	snp_out->var_id = (char**) malloc(length * sizeof(char*));
	snp_out->chr = (char**) malloc(length * sizeof(char*));
	snp_out->pos_morgans = (double*) malloc(length * sizeof(double));
	snp_out->pos = (uint64_t*) malloc(length * sizeof(uint64_t));
	snp_out->ref = (char**) malloc(length * sizeof(char*));
	snp_out->alt = (char**) malloc(length * sizeof(char*));
	snp_out->rev_idx = (id_map_str*) malloc(sizeof(id_map_str));
	snp_out->hash = 0;
	snp_out->rev_idx->map = kh_init(ID_MAP_STR);
	size_t i = 0;
	STAILQ_FOREACH(tmp_node, head, nodes) {
		// assign values
		snp_out->var_id[i] = strdup(snp_in->var_id[tmp_node->idx]);
		snp_out->chr[i] = strdup(snp_in->chr[tmp_node->idx]);
		snp_out->pos_morgans[i] = snp_in->pos_morgans[tmp_node->idx];
		snp_out->pos[i] = snp_in->pos[tmp_node->idx];
		snp_out->ref[i] = strdup(snp_in->ref[tmp_node->idx]);
		snp_out->alt[i] = strdup(snp_in->alt[tmp_node->idx]);
		// calculate hash
		snp_out->hash *= 17;
		snp_out->hash = snp_out->hash ^ hash_str(snp_out->var_id[i]);
		// reverse index hash table management
		int ret;
		khiter_t k = kh_put(ID_MAP_STR, snp_out->rev_idx->map, snp_out->var_id[i], &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert variant name into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Multiple instances of variant %s found.\n", snp_out->var_id[i]);
			exit(EXIT_FAILURE);
		}
		kh_value(snp_out->rev_idx->map, k) = i;
		i++;
	}
	return 0;
}

short filter_ind_data(ind_data* ind_in, ind_data* ind_out, struct idx_head* head) {
	// get length
	size_t length = 0;
	struct idx_node* tmp_node;
	STAILQ_FOREACH(tmp_node, head, nodes) {
		length++;
	}
	if(length == 0) {
		return -1;
	}
	// make struct
	ind_out->length = length,
	ind_out->ind_id = (char**) malloc(length * sizeof(char*));
	ind_out->sex = (char**) malloc(length * sizeof(char*));
	ind_out->population = (char**) malloc(length * sizeof(char*));
	ind_out->rev_idx = (id_map_ind*) malloc(sizeof(id_map_ind));
	ind_out->rev_idx->map = kh_init(ID_MAP_IND);	
	size_t i = 0;
	STAILQ_FOREACH(tmp_node, head, nodes) {
		// assign values
		ind_out->ind_id[i] = strdup(ind_in->ind_id[tmp_node->idx]);
		ind_out->sex[i] = strdup(ind_in->sex[tmp_node->idx]);
		ind_out->population[i] = strdup(ind_in->population[tmp_node->idx]);
		// calculate hash
		ind_out->hash *= 17;
		ind_out->hash = ind_out->hash ^ hash_str(ind_out->ind_id[i]);
		// reverse index hash table management
		int ret;
		ind_idx ind_idx_struct = (ind_idx){ind_out->ind_id[i], ind_out->population[i]};
		khiter_t k = kh_put(ID_MAP_IND, ind_out->rev_idx->map, ind_idx_struct, &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert individual into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Multiple instances of individual (%s, %s) found.\n", ind_out->ind_id[i], ind_out->population[i]);
			exit(EXIT_FAILURE);
		}
		kh_value(ind_out->rev_idx->map, k) = i;
		i++;
	}
	return 0;
}

size_t intersect_snp_data(snp_data* snp1, snp_data* snp2, struct idx_head* head1, struct idx_head* head2) {
	size_t length = 0;
	for(size_t idx1 = 0; idx1 < snp1->length; idx1++) {
		size_t idx2;
		short ret = get_snp_idx(snp2, snp1->var_id[idx1], &idx2);
		if(ret == 0) {
			struct idx_node* idn1 = (struct idx_node*)malloc(sizeof(struct idx_node));
			idn1->idx = idx1;
			STAILQ_INSERT_TAIL(head1, idn1, nodes);
			struct idx_node* idn2 = (struct idx_node*)malloc(sizeof(struct idx_node));
			idn1->idx = idx2;
			STAILQ_INSERT_TAIL(head2, idn2, nodes);
			length++;
		}
	}
	return length;
}

void append_ind_data(ind_data* ind1, ind_data* ind2, ind_data* ind_out) {
	ind_out->length = ind1->length + ind2->length;
	ind_out->hash = 0;
	ind_out->ind_id = (char**) malloc(ind_out->length * sizeof(char*));
	ind_out->sex = (char**) malloc(ind_out->length * sizeof(char*));
	ind_out->population = (char**) malloc(ind_out->length * sizeof(char*));
	// copy ind1 data
	memcpy(ind_out->ind_id, ind1->ind_id, ind1->length);
	memcpy(ind_out->sex, ind1->sex, ind1->length);
	memcpy(ind_out->population, ind1->population, ind1->length);
	// copy ind2 data
	memcpy(ind_out->ind_id + ind1->length, ind2->ind_id, ind2->length);
	memcpy(ind_out->sex + ind1->length, ind2->sex, ind2->length);
	memcpy(ind_out->population + ind1->length, ind2->population, ind2->length);
	// init hash
	ind_out->rev_idx = (id_map_ind*) malloc(sizeof(id_map_ind));
	ind_out->rev_idx->map = kh_init(ID_MAP_IND);
	for(size_t idx = 0; idx < ind_out->length; idx++) {
		int ret;
		ind_idx ind_idx_struct = (ind_idx){ind_out->ind_id[idx], ind_out->population[idx]};
		khiter_t k = kh_put(ID_MAP_IND, ind_out->rev_idx->map, ind_idx_struct, &ret);
		if(ret == -1) {
			fprintf(stderr, "Error: Failed to insert individual into hash table!\n");
			exit(EXIT_FAILURE);
		} else if(ret == 0) {
			fprintf(stderr, "Error: Instances of individual (%s, %s) found in both individual objects.\n", ind_out->ind_id[idx], ind_out->population[idx]);
			exit(EXIT_FAILURE);
		}
		kh_value(ind_out->rev_idx->map, k) = idx;
		ind_out->hash *= 17;
		ind_out->hash = ind_out->hash ^ hash_str(ind_out->ind_id[idx]);
	}
}

pam_file_reader pam_file_reader_init(char* filename, snp_data* snp_info, ind_data* ind_info) {
	pam_file_reader pf;
	size_t file_size = get_filesize(filename);
	pf.is_hdr_read = false;
	pf.is_open = true;
	pf.idx = 0;
	pf.n_ind = ind_info->length;
	pf.n_snp = snp_info->length;
	pf.record_size = file_size / (snp_info->length + 1);
	if ((pf.record_size * (snp_info->length + 1)) != file_size) {
		fprintf(stderr, "Invalid PACKEDANCESTRYMAP file. File size must be a multiple of the number of SNPs.\n");
		exit(EXIT_FAILURE);
	}
	pf.fp = safe_read(filename, "rb");
	return pf;
}

egn_file_reader egn_file_reader_init(char* filename, snp_data* snp_info, ind_data* ind_info) {
	egn_file_reader ef;
	ef.is_open = true;
	ef.n_ind = ind_info->length;
	ef.n_snp = snp_info->length;
	size_t file_size = get_filesize(filename);
	ef.idx = 0;
	if (((ind_info->length + 1) * snp_info->length) != file_size) {
		fprintf(stderr, "ERROR: Invalid EIGENSTRAT file. Expected file of size %zu bytes.\n", (ind_info->length + 1) * snp_info->length);
		exit(EXIT_FAILURE);
	}
	ef.fp = safe_read(filename, "rb");
	return ef;
}

int close_pam_file_reader(pam_file_reader* pf) {
	pf->is_open = false;
	return fclose(pf->fp);
}

int close_egn_file_reader(egn_file_reader* ef) {
	ef->is_open = false;
	return fclose(ef->fp);
}

int close_pam_file_writer(pam_file_writer* pfw) {
	pfw->is_open = false;
	return fclose(pfw->fp);
}

int close_egn_file_writer(egn_file_writer* efw) {
	efw->is_open = false;
	return fclose(efw->fp);
}

hdr_data read_pam_header(pam_file_reader* pf) {
	if(pf->is_hdr_read) {
		fprintf(stderr, "Header already read!");
		exit(EXIT_FAILURE);
	}
	if(!pf->is_open) {
		fprintf(stderr, "ERROR: pam_file_reader has been closed.\n");
		exit(EXIT_FAILURE);
	}
	hdr_data hdr_info;
	char* hdr = (char*) malloc(sizeof(char) * pf->record_size);
	fread(hdr, 1, pf->record_size, pf->fp);
	sscanf(hdr, "GENO   %zu %zu %x %x", &hdr_info.n_ind, &hdr_info.n_snp, &hdr_info.ind_hash, &hdr_info.snp_hash);
	free(hdr);
	pf->is_hdr_read = true;
	return hdr_info;
}

uint8_t* read_pam_record(pam_file_reader* pf) {
	if(!pf->is_hdr_read) {
		fprintf(stderr, "Header must be read before reading records!\n");
		exit(EXIT_FAILURE);
	}
	if(!pf->is_open) {
		fprintf(stderr, "ERROR: pam_file_reader has been closed.\n");
		exit(EXIT_FAILURE);
	}
	if(pf->idx == pf->n_snp) {
		return NULL;
	}
	uint8_t* record = (uint8_t*)malloc(pf->n_ind * sizeof(uint8_t));
	size_t num_leftover_bytes = pf->record_size - (int)ceil((float)(pf->n_ind*RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE);
	uint8_t record_byte;
	for(size_t i = 0; i < pf->n_ind; i++) {
		uint8_t elem_pos = i%RECORD_ELEMS_PER_BYTE;
		uint8_t shift_by = (BITS_IN_BYTE-RECORD_ELEM_SIZE_BITS) - (RECORD_ELEM_SIZE_BITS*elem_pos);
		if(elem_pos == 0) {
			record_byte = getc(pf->fp);
		}
		record[i] = (record_byte & (RECORD_ELEMS_MASK_BASE << shift_by)) >> shift_by;
	}
	pf->idx += 1;
	fseek(pf->fp, num_leftover_bytes, SEEK_CUR);  // skip useless bytes
	return record;
}

uint8_t* read_egn_record(egn_file_reader* ef) {
	if(!ef->is_open) {
		fprintf(stderr, "ERROR: egn_file_reader has been closed.\n");
		exit(EXIT_FAILURE);
	}
	if(ef->idx == ef->n_snp) {
		return NULL;
	}
	uint8_t* record = (uint8_t*)malloc(ef->n_ind * sizeof(uint8_t));
	char* line = NULL;
	size_t size = 0;
	ssize_t nread = getline(&line, &size, ef->fp);
	if (nread != (ef->n_ind+1)) {
		fprintf(stderr, "ERROR: improperly formatted EIGENSTRAT geno file. Expected %zu entries in line, got %li.\n", ef->n_ind+1, nread);
		exit(EXIT_FAILURE);
	}
	for(size_t i = 0; i < ef->n_ind; i++) {
		// check if valid character
		if(!((line[i] == '0') || (line[i] == '1') || (line[i] == '2') || (line[i] == '9'))) {
			fprintf(stderr, "ERROR: characters must be '0', '1', '2', or '9'.\n");
			exit(EXIT_FAILURE);
		}
		record[i] = ((line[i] == '0') * 0) + ((line[i] == '1') * 1) + ((line[i] == '2') * 2) + ((line[i] == '9') * NAN_VAL);
	}
	free(line);
	ef->idx += 1;
	return record;
}

short goto_var_egn(egn_file_reader* ef, snp_data* snp_info, char* var_name) {
	if(!ef->is_open) {
		fprintf(stderr, "ERROR: egn_file_reader closed.\n");
		exit(EXIT_FAILURE);
	}
	size_t idx_go;
	short ret = get_snp_idx(snp_info, var_name, &idx_go);
	if(ret == 0) {
		fseek(ef->fp, idx_go * (ef->n_ind+1), SEEK_SET);
		ef->idx = idx_go;
	}
	return ret;
}

short goto_var_pam(pam_file_reader* pf, snp_data* snp_info, char* var_name) {
	if(!pf->is_open) {
		fprintf(stderr, "ERROR: pam_file_reader closed.\n");
		exit(EXIT_FAILURE);
	}
	if(!pf->is_hdr_read) {
		fprintf(stderr, "ERROR: header must be read before going to variant.\n");
		exit(EXIT_FAILURE);	
	}
	size_t idx_go;
	short ret = get_snp_idx(snp_info, var_name, &idx_go);
	if(ret == 0) {
		fseek(pf->fp, (idx_go + 1) * (pf->record_size), SEEK_SET);
		pf->idx = idx_go;
	}
	return ret;
}

pam_file_writer pam_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info) {
	pam_file_writer pfw;
	pfw.hdr_written = false;
	pfw.is_open = true;
	pfw.n_snp = snp_info->length;
	pfw.n_ind = ind_info->length;
	pfw.n_written_snp = 0;
	pfw.fp = safe_read(filename, "wb+");
	return pfw;
}

egn_file_writer egn_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info) {
	egn_file_writer efw;
	efw.is_open = true;
	efw.n_snp = snp_info->length;
	efw.n_ind = ind_info->length;
	efw.n_written_snp = 0;
	efw.fp = safe_read(filename, "w+");
	return efw;
}

void write_pam_header(pam_file_writer* pfw, snp_data* snp_info, ind_data* ind_info) {
	if(pfw->hdr_written) {
		fprintf(stderr, "ERROR: pam_file_writer object header has already been written!\n");
		exit(EXIT_FAILURE);
	}
	if(!pfw->is_open) {
		fprintf(stderr, "ERROR: PACKEDANCESTRYMAP file is closed!\n");
		exit(EXIT_FAILURE);
	}
	int num_chars = fprintf(pfw->fp, "GENO   %zu %zu %x %x", 	pfw->n_ind, pfw->n_snp, ind_info->hash, snp_info->hash);
	if(num_chars < 0) {
		fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
		exit(EXIT_FAILURE);
	}
	size_t record_size = MAX(num_chars, (int)ceil((float)(pfw->n_ind * RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE));
  	size_t n_trailing_bytes_hdr	= record_size - num_chars;
	for(size_t i = 0; i < n_trailing_bytes_hdr; i++) {
		int ret = fputc('\0', pfw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
			exit(EXIT_FAILURE);
		}
	}
	pfw->hdr_written = true;
	pfw->record_size = record_size;
}

void write_pam_record(pam_file_writer* pfw, uint8_t* dosages) {
	if(!pfw->hdr_written) {
		fprintf(stderr, "ERROR: header has not been written for PACKEDANCESTRYMAP file!\n");
		exit(EXIT_FAILURE);
	}
	if(!pfw->is_open) {
		fprintf(stderr, "ERROR: PACKEDANCESTRYMAP file is closed!\n");
		exit(EXIT_FAILURE);
	}
	if(pfw->n_written_snp >= pfw->n_snp) {
		fprintf(stderr, "ERROR: Cannot write to PACKEDANCESTRYMAP file. All SNP records have been written.\n");
		exit(EXIT_FAILURE);
	}
	uint8_t record_byte = 0;
	for(size_t i = 0; i < pfw->n_ind; i++) {
		if((i != 0) && ((i%(BITS_IN_BYTE/RECORD_ELEM_SIZE_BITS)) == 0)) {
			int ret = fputc(record_byte, pfw->fp);
			if(ret == EOF) {
				fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
				exit(EXIT_FAILURE);
			}
			record_byte = 0;
		}
		uint8_t elem_pos = i%RECORD_ELEMS_PER_BYTE;
		uint8_t shift_by = (BITS_IN_BYTE-RECORD_ELEM_SIZE_BITS) - (RECORD_ELEM_SIZE_BITS*elem_pos);
		record_byte = record_byte | ((RECORD_ELEMS_MASK_BASE << shift_by) & (dosages[i] << shift_by));
	}

	// write last record if not written
	if((pfw->n_ind % (BITS_IN_BYTE/RECORD_ELEM_SIZE_BITS)) != 0) {
		int ret = fputc(record_byte, pfw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
			exit(EXIT_FAILURE);
		}
	}

	size_t n_trailing_bytes = pfw->record_size - (int)ceil((float)(pfw->n_ind * RECORD_ELEM_SIZE_BITS) / BITS_IN_BYTE);
	// write trailing bytes
	for(size_t i = 0; i < n_trailing_bytes; i++) {
		int ret = fputc('\0', pfw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to PACKEDANCESTRYMAP file.\n");
			exit(EXIT_FAILURE);
		}
	}
	pfw->n_written_snp += 1;
}

void write_egn_record(egn_file_writer* efw, uint8_t* dosages) {
	if(!efw->is_open) {
		fprintf(stderr, "ERROR: EIGENSTRAT file is closed!\n");
		exit(EXIT_FAILURE);
	}
	if(efw->n_written_snp >= efw->n_snp) {
		fprintf(stderr, "ERROR: Cannot write to EIGENSTRAT file. All SNP records have been written.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i = 0; i < efw->n_ind; i++) {
		char dsg = ('0'*(dosages[i] == 0)) + ('1'*(dosages[i] == 1)) + ('2'*(dosages[i] == 2)) + ('9'*(dosages[i] == NAN_VAL));
		int ret = fputc(dsg, efw->fp);
		if(ret == EOF) {
			fprintf(stderr, "ERROR: unsuccessful write to EIGENSTRAT file.\n");
			exit(EXIT_FAILURE);
		}
	}
	int ret = fputc('\n', efw->fp);
	if(ret == EOF) {
		fprintf(stderr, "ERROR: unsuccessful write to EIGENSTRAT file.\n");
		exit(EXIT_FAILURE);
	}
	efw->n_written_snp += 1;
}

void free_str_array(char** arr, size_t length) {
	for(size_t i = 0; i < length; i++) {
		free(arr[i]);
	}
	free(arr);
}

void free_snp_data(snp_data* snp_info) {
	free_str_array(snp_info->var_id, snp_info->length);
	free_str_array(snp_info->chr, snp_info->length);
	free(snp_info->pos_morgans);
	free(snp_info->pos);
	free_str_array(snp_info->ref, snp_info->length);
	free_str_array(snp_info->alt, snp_info->length);
	kh_destroy(ID_MAP_STR, snp_info->rev_idx->map);
	free(snp_info->rev_idx);
}

void free_ind_data(ind_data* ind_info) {
	free_str_array(ind_info->ind_id, ind_info->length);
	free_str_array(ind_info->sex, ind_info->length);
	free_str_array(ind_info->population, ind_info->length);
	kh_destroy(ID_MAP_IND, ind_info->rev_idx->map);
	free(ind_info->rev_idx);
}

void free_idx_list(struct idx_head* head) {
	while(!STAILQ_EMPTY(head)) {
		struct idx_node* tmp = STAILQ_FIRST(head);
		STAILQ_REMOVE_HEAD(head, nodes);
		free(tmp);
	}
}

void free_str_list(struct str_list_head* head) {
	while(!STAILQ_EMPTY(head)) {
		struct str_node* tmp = STAILQ_FIRST(head);
		STAILQ_REMOVE_HEAD(head, nodes);
		free(tmp->str);
		free(tmp);
	}
}

void free_ind_idx_list(struct ind_idx_head* head) {
	while(!STAILQ_EMPTY(head)) {
		struct ind_idx_node* tmp = STAILQ_FIRST(head);
		STAILQ_REMOVE_HEAD(head, nodes);
		free(tmp->iidx->ind_id);
		free(tmp->iidx->ind_pop);
		free(tmp->iidx);
		free(tmp);
	}
}
