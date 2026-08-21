#define  _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/queue.h>
#include "unity.h"
#include "dotgeno.h"

void setUp(void) { }

void tearDown(void) {
	remove("/tmp/tmppam.geno");
	remove("/tmp/tmpegn.geno");
	remove("/tmp/tmptgn.tgeno");

}

size_t get_random_index(size_t length) {
	srand(time(NULL));
	return((size_t)rand() % length);
}

void test_equality_egn_pam(void) {
	char* file_base_egn = "test/data/PAM_EGN/egn_example_%i.%s";
	char* file_base_pam = "test/data/PAM_EGN/pam_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base_pam, i, "snp");
		snp_data snp_info_pam = read_snp_file(buf);
		sprintf(buf, file_base_pam, i, "ind");
		ind_data ind_info_pam = read_ind_file(buf);
		sprintf(buf, file_base_pam, i, "geno");
		pam_file_reader pfr = pam_file_reader_init(buf, &snp_info_pam, &ind_info_pam);
		read_pam_header(&pfr);

		sprintf(buf, file_base_egn, i, "snp");
		snp_data snp_info_egn = read_snp_file(buf);
		sprintf(buf, file_base_egn, i, "ind");
		ind_data ind_info_egn = read_ind_file(buf);
		sprintf(buf, file_base_egn, i, "geno");
		egn_file_reader efr = egn_file_reader_init(buf, &snp_info_egn, &ind_info_egn);

		TEST_ASSERT_EQUAL_UINT(snp_info_pam.length, snp_info_egn.length);
		TEST_ASSERT_EQUAL_UINT(ind_info_pam.length, ind_info_egn.length);
		
		uint8_t* pam_record;
		uint8_t* egn_record;
		while(1) {
			pam_record = read_pam_record(&pfr);
			egn_record = read_egn_record(&efr);
			if(pam_record == NULL) {
				if(egn_record) {
					free(egn_record);
					TEST_FAIL();
				}
				break;
			}
			TEST_ASSERT_EQUAL_UINT8_ARRAY(pam_record, egn_record, ind_info_pam.length);
			free(pam_record);
			free(egn_record);
		}
		free_snp_data(&snp_info_pam);
		free_snp_data(&snp_info_egn);
		free_ind_data(&ind_info_pam);
		free_ind_data(&ind_info_egn);
		close_pam_file_reader(&pfr);
		close_egn_file_reader(&efr);
	}
}

void test_goto_equality_pam_egn(void) {
	char* file_base_egn = "test/data/PAM_EGN/egn_example_%i.%s";
	char* file_base_pam = "test/data/PAM_EGN/pam_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base_pam, i, "snp");
		snp_data snp_info_pam = read_snp_file(buf);
		sprintf(buf, file_base_pam, i, "ind");
		ind_data ind_info_pam = read_ind_file(buf);
		sprintf(buf, file_base_pam, i, "geno");
		pam_file_reader pfr = pam_file_reader_init(buf, &snp_info_pam, &ind_info_pam);
		read_pam_header(&pfr);

		sprintf(buf, file_base_egn, i, "snp");
		snp_data snp_info_egn = read_snp_file(buf);
		sprintf(buf, file_base_egn, i, "ind");
		ind_data ind_info_egn = read_ind_file(buf);
		sprintf(buf, file_base_egn, i, "geno");
		egn_file_reader efr = egn_file_reader_init(buf, &snp_info_egn, &ind_info_egn);

		TEST_ASSERT_EQUAL_UINT(snp_info_pam.length, snp_info_egn.length);
		TEST_ASSERT_EQUAL_UINT(ind_info_pam.length, ind_info_egn.length);
		
		int n_iters = 7;
		uint8_t* pam_record;
		uint8_t* egn_record;
		for(int j = 0; j < n_iters; j++) {
			short ret;
			ret = goto_var_pam(&pfr, &snp_info_pam, snp_info_pam.var_id[get_random_index(snp_info_pam.length)]);
			if(ret == -1) { TEST_FAIL(); };
			ret = goto_var_egn(&efr, &snp_info_egn, snp_info_egn.var_id[get_random_index(snp_info_egn.length)]);
			if(ret == -1) { TEST_FAIL(); };
			pam_record = read_pam_record(&pfr);
			egn_record = read_egn_record(&efr);
			TEST_ASSERT_EQUAL_UINT8_ARRAY(pam_record, egn_record, ind_info_pam.length);
			free(pam_record);
			free(egn_record);
		}
		free_snp_data(&snp_info_pam);
		free_snp_data(&snp_info_egn);
		free_ind_data(&ind_info_pam);
		free_ind_data(&ind_info_egn);
		close_pam_file_reader(&pfr);
		close_egn_file_reader(&efr);
	}
}

void test_ind_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/total/ind", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char ind_file[37];
		uint32_t exp_hash;
		sscanf(buf, "%s %x\n", ind_file, &exp_hash);
		ind_data ind_info = read_ind_file(ind_file);
		TEST_ASSERT_EQUAL_UINT32(exp_hash, ind_info.hash);
		free_ind_data(&ind_info);
	}
	free(buf);
}

void test_snp_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/total/snp", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		uint32_t exp_hash;
		sscanf(buf, "%s %x\n", snp_file, &exp_hash);
		snp_data snp_info = read_snp_file(snp_file);
		TEST_ASSERT_EQUAL_UINT32(exp_hash, snp_info.hash);
		free_snp_data(&snp_info);
	}
	free(buf);
}

void test_ind_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/ind", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char ind_file[37];
		uint32_t exp_hash;
		char* idx_str = malloc(sizeof(char) * ((bufsize - 46) + 1));
		sscanf(buf, "%s\t%s\t%x\n", ind_file, idx_str, &exp_hash);
		ind_data ind_info = read_ind_file(ind_file);
		
		// init idx list
		struct idx_head head;
		TAILQ_INIT(&head);
		char* elem = strtok(idx_str, ","); 
		// add elems
		while(elem) {
			struct idx_node* node = malloc(sizeof(struct idx_node));
			node->idx = atoi(elem);
			TAILQ_INSERT_TAIL(&head, node, nodes);
			elem = strtok(NULL, ",");
		}
		ind_data ind_filt;
		short ret = filter_ind_data(&ind_info, &ind_filt, &head);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, ind_filt.hash);
		free_ind_data(&ind_info);
		free(idx_str);
		free_idx_list(&head);
	}
	free(buf);
}

void test_snp_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/snp", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		uint32_t exp_hash;
		char* idx_str = malloc(sizeof(char) * ((bufsize - 46) + 1));
		sscanf(buf, "%s\t%s\t%x\n", snp_file, idx_str, &exp_hash);
		snp_data snp_info = read_snp_file(snp_file);
		
		// init idx list
		struct idx_head head;
		TAILQ_INIT(&head);
		char* elem = strtok(idx_str, ","); 
		// add elems
		while(elem) {
			struct idx_node* node = malloc(sizeof(struct idx_node));

			node->idx = atoi(elem);
			TAILQ_INSERT_TAIL(&head, node, nodes);
			elem = strtok(NULL, ",");
		}
		snp_data snp_filt;
		short ret = filter_snp_data(&snp_info, &snp_filt, &head);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, snp_filt.hash);
		free_snp_data(&snp_info);
		free(idx_str);
		free_idx_list(&head);
	}
	free(buf);
}

void test_pop_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/pop", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char ind_file[37];
		uint32_t exp_hash;
		char pop_list[1000];
		sscanf(buf, "%s\t%s\t%x\n", ind_file, pop_list, &exp_hash);
		ind_data ind_info = read_ind_file(ind_file);
		
		// init pop array
		// get number of elements
		size_t n_elems = 1;
		int len = strlen(pop_list);
		for(int i = 0; i < len; i++) {
			if(pop_list[i] == ',') { n_elems++; }
		}
		char** pops = (char**)malloc(n_elems * sizeof(char*));
		char* token = strtok(pop_list, ",");
		pops[0] = token;
		size_t cur_i = 1;
		while (token != NULL) {
			token = strtok(NULL, ",");
			if(token) {
				pops[cur_i] = token;
				cur_i++;
			}	
		}

		struct idx_head head_idx;
		TAILQ_INIT(&head_idx);
		get_multiple_pops(&ind_info, pops, n_elems, &head_idx, NULL);

		ind_data ind_filt;
		short ret = filter_ind_data(&ind_info, &ind_filt, &head_idx);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, ind_filt.hash);
		free_ind_data(&ind_info);
		free_ind_data(&ind_filt);
		free(pops);
		free_idx_list(&head_idx);
	}
	free(buf);
}

void test_chr_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/chr", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		uint32_t exp_hash;
		char chr_list[1000];
		sscanf(buf, "%s\t%s\t%x\n", snp_file, chr_list, &exp_hash);
		snp_data snp_info = read_snp_file(snp_file);
		
		// init chr array
		// get number of elements
		size_t n_elems = 1;
		int len = strlen(chr_list);
		for(int i = 0; i < len; i++) {
			if(chr_list[i] == ',') { n_elems++; }
		}
		char** chrs = (char**)malloc(n_elems * sizeof(char*));
		char* token = strtok(chr_list, ",");
		chrs[0] = token;
		size_t cur_i = 1;
		while (token != NULL) {
			token = strtok(NULL, ",");
			if(token) {
				chrs[cur_i] = token;
				cur_i++;
			}	
		}
		struct idx_head head_idx;
		TAILQ_INIT(&head_idx);
		get_multiple_chrs(&snp_info, chrs, n_elems, &head_idx);

		snp_data snp_filt;
		short ret = filter_snp_data(&snp_info, &snp_filt, &head_idx);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, snp_filt.hash);
		free_snp_data(&snp_info);
		free_snp_data(&snp_filt);
		free(chrs);
		free_idx_list(&head_idx);
	}
	free(buf);
}

void test_range_filtered_hash(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/hash_test/filtered/range", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		char chr[3];
		uint32_t exp_hash;
		uint64_t start;
		uint64_t end;
		sscanf(buf, "%s\t%s\t%zu\t%zu\t%x\n", snp_file, chr, &start, &end, &exp_hash);
		snp_data snp_info = read_snp_file(snp_file);
		char* chrs[1];
		uint64_t starts[1];
		uint64_t ends[1];
		chrs[0] = chr;
		starts[0] = start;
		ends[0] = end;

		struct idx_head head_idx;
		TAILQ_INIT(&head_idx);
		get_multiple_ranges(&snp_info, chrs, starts, ends, 1, &head_idx);

		snp_data snp_filt;
		short ret = filter_snp_data(&snp_info, &snp_filt, &head_idx);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_UINT32(exp_hash, snp_filt.hash);
		free_snp_data(&snp_info);
		free_snp_data(&snp_filt);
		free_idx_list(&head_idx);
	}
	free(buf);
}

void test_read_write_eqv_pam(void) {
	char* file_base = "test/data/PAM_EGN/pam_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base, i, "snp");
		snp_data snp_info = read_snp_file(buf);
		sprintf(buf, file_base, i, "ind");
		ind_data ind_info = read_ind_file(buf);
		sprintf(buf, file_base, i, "geno");
		pam_file_reader pfr = pam_file_reader_init(buf, &snp_info, &ind_info);
		read_pam_header(&pfr);

		pam_file_writer pfw = pam_file_writer_init("/tmp/tmppam.geno", &snp_info, &ind_info);
		write_pam_header(&pfw, &snp_info, &ind_info);

		uint8_t* record;
		while(record = read_pam_record(&pfr)) {
			write_pam_record(&pfw, record);
			free(record);
		}
		TEST_ASSERT_EQUAL_INT(pfr.n_snp, pfw.n_written_snp);
		close_pam_file_writer(&pfw);
		close_pam_file_reader(&pfr);

		pam_file_reader pfr_in = pam_file_reader_init(buf, &snp_info, &ind_info);
		pam_file_reader pfr_out = pam_file_reader_init("/tmp/tmppam.geno", &snp_info, &ind_info);
		TEST_ASSERT_EQUAL_INT(pfr_in.n_snp, pfr_out.n_snp);
		TEST_ASSERT_EQUAL_INT(pfr_in.n_ind, pfr_out.n_ind);
		
		hdr_data hdr1 = read_pam_header(&pfr_in);
		hdr_data hdr2 = read_pam_header(&pfr_out);
		TEST_ASSERT_EQUAL_INT(hdr1.n_ind, hdr2.n_ind);
		TEST_ASSERT_EQUAL_INT(hdr1.n_snp, hdr2.n_snp);
		TEST_ASSERT_EQUAL_UINT32(hdr1.ind_hash, hdr2.ind_hash);
		TEST_ASSERT_EQUAL_UINT32(hdr1.snp_hash, hdr2.snp_hash);
		uint8_t* record1;
		uint8_t* record2;
		while(1) {
			record1 = read_pam_record(&pfr_in);
			record2 = read_pam_record(&pfr_out);
			if(record1 == NULL) {
				if(record2) {
					TEST_FAIL();
				}
				break;
			}
			TEST_ASSERT_EQUAL_UINT8_ARRAY(record1, record2, ind_info.length);
			free(record1);
			free(record2);
		}
		free_snp_data(&snp_info);
		free_ind_data(&ind_info);
	}
}

void test_read_write_eqv_tgn(void) {
	char* file_base = "test/data/TGN/tgn_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base, i, "snp");
		snp_data snp_info = read_snp_file(buf);
		sprintf(buf, file_base, i, "ind");
		ind_data ind_info = read_ind_file(buf);
		sprintf(buf, file_base, i, "tgeno");
		tgn_file_reader tfr = tgn_file_reader_init(buf, &snp_info, &ind_info);
		read_tgn_header(&tfr);

		tgn_file_writer tfw = tgn_file_writer_init("/tmp/tmptgn.tgeno", &snp_info, &ind_info);
		write_tgn_header(&tfw, &snp_info, &ind_info);

		uint8_t* record;
		while(record = read_tgn_record(&tfr)) {
			write_tgn_record(&tfw, record);
			free(record);
		}
		TEST_ASSERT_EQUAL_INT(tfr.n_ind, tfw.n_written_ind);
		close_tgn_file_writer(&tfw);
		close_tgn_file_reader(&tfr);

		tgn_file_reader tfr_in = tgn_file_reader_init(buf, &snp_info, &ind_info);
		tgn_file_reader tfr_out = tgn_file_reader_init("/tmp/tmptgn.tgeno", &snp_info, &ind_info);
		TEST_ASSERT_EQUAL_INT(tfr_in.n_snp, tfr_out.n_snp);
		TEST_ASSERT_EQUAL_INT(tfr_in.n_ind, tfr_out.n_ind);
		
		hdr_data hdr1 = read_tgn_header(&tfr_in);
		hdr_data hdr2 = read_tgn_header(&tfr_out);
		TEST_ASSERT_EQUAL_INT(hdr1.n_ind, hdr2.n_ind);
		TEST_ASSERT_EQUAL_INT(hdr1.n_snp, hdr2.n_snp);
		TEST_ASSERT_EQUAL_UINT32(hdr1.ind_hash, hdr2.ind_hash);
		TEST_ASSERT_EQUAL_UINT32(hdr1.snp_hash, hdr2.snp_hash);
		uint8_t* record1;
		uint8_t* record2;
		while(1) {
			record1 = read_tgn_record(&tfr_in);
			record2 = read_tgn_record(&tfr_out);
			if(record1 == NULL) {
				if(record2) {
					TEST_FAIL();
				}
				break;
			}
			TEST_ASSERT_EQUAL_UINT8_ARRAY(record1, record2, snp_info.length);
			free(record1);
			free(record2);
		}
		free_snp_data(&snp_info);
		free_ind_data(&ind_info);
	}
}

void test_read_write_eqv_egn(void) {
	char* file_base = "test/data/PAM_EGN/egn_example_%i.%s";
	char buf[38];
    int n_sets = 20;
	for(int i = 1; i <= n_sets; i++) {
		sprintf(buf, file_base, i, "snp");
		snp_data snp_info = read_snp_file(buf);
		sprintf(buf, file_base, i, "ind");
		ind_data ind_info = read_ind_file(buf);
		sprintf(buf, file_base, i, "geno");
		egn_file_reader pfr = egn_file_reader_init(buf, &snp_info, &ind_info);

		egn_file_writer pfw = egn_file_writer_init("/tmp/tmpegn.geno", &snp_info, &ind_info);

		uint8_t* record;
		while(record = read_egn_record(&pfr)) {
			write_egn_record(&pfw, record);
			free(record);
		}
		TEST_ASSERT_EQUAL_INT(pfr.n_snp, pfw.n_written_snp);
		close_egn_file_writer(&pfw);
		close_egn_file_reader(&pfr);

		egn_file_reader pfr_in = egn_file_reader_init(buf, &snp_info, &ind_info);
		egn_file_reader pfr_out = egn_file_reader_init("/tmp/tmpegn.geno", &snp_info, &ind_info);
		TEST_ASSERT_EQUAL_INT(pfr_in.n_snp, pfr_out.n_snp);
		TEST_ASSERT_EQUAL_INT(pfr_in.n_ind, pfr_out.n_ind);
		
		uint8_t* record1;
		uint8_t* record2;
		while(1) {
			record1 = read_egn_record(&pfr_in);
			record2 = read_egn_record(&pfr_out);
			if(record1 == NULL) {
				if(record2) {
					TEST_FAIL();
				}
				break;
			}
			TEST_ASSERT_EQUAL_UINT8_ARRAY(record1, record2, ind_info.length);
			free(record1);
			free(record2);
		}
		free_snp_data(&snp_info);
		free_ind_data(&ind_info);
	}
}

void test_get_idx_snp(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/single_index_test/snp", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		char var_name[100];
		size_t idx_exp;
		sscanf(buf, "%s\t%s\t%lu\n", snp_file, var_name, &idx_exp);
		snp_data snp_info = read_snp_file(snp_file);
		size_t idx_actual;
		short ret = get_snp_idx(&snp_info, var_name, &idx_actual);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_INT(idx_exp, idx_actual);
		// now test that a gibberish var_name fails
		ret = get_snp_idx(&snp_info, "GIBBERISHGIBBERGABER", &idx_actual);
		if(ret == 0) { TEST_FAIL(); }
	}
	free(buf);
}

void test_get_idx_ind(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/single_index_test/ind", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char ind_file[37];
		char idx_str[100];
		size_t idx_exp;
		sscanf(buf, "%s\t%s\t%lu\n", ind_file, idx_str, &idx_exp);

		// split string
		int len_idx_str = strlen(idx_str);
		int split_pos;
		for(int i = 0; i < len_idx_str; i++) {
			if(idx_str[i] == ',') {
				split_pos = i;
				idx_str[i] = '\0';
			}
		}
		// store info here
		char* ind_id = idx_str;
		char* ind_pop = &idx_str[split_pos + 1];	
		ind_data ind_info = read_ind_file(ind_file);
		size_t idx_actual;
		short ret = get_ind_idx(&ind_info, ind_id, ind_pop, &idx_actual);
		if(ret == -1) { TEST_FAIL(); }
		TEST_ASSERT_EQUAL_INT(idx_exp, idx_actual);
		// now test that a gibberish var_name fails
		ret = get_ind_idx(&ind_info, "GUhjdshjdchjdhj", "hjfdhjjhj", &idx_actual);
		if(ret == 0) { TEST_FAIL(); }
	}
	free(buf);
}

void test_multiple_index_snp(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/multiple_index_test/snp", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char snp_file[37];
		size_t n_elems;
		sscanf(buf, "%s\t%lu\t%*s\t%*s\n", snp_file, &n_elems);
		char** snp_names = malloc(sizeof(char*) * n_elems);
		size_t* exp_idx = malloc(sizeof(size_t) * n_elems);
		char snp_str[10000];
		char idx_str[10000];
		sscanf(buf, "%*s\t%*lu\t%s\t%s\n", snp_str, idx_str);
		
		char* elem = strtok(snp_str, ",");
		size_t i = 0;
		while(elem) {
			snp_names[i] = strdup(elem);
			elem = strtok(NULL, ",");
			i++;
		}
		if(strcmp(idx_str, "NA") != 0) {
			elem = strtok(idx_str, ",");
			i = 0;
			while(elem) {
				exp_idx[i] = atoi(elem);
				elem = strtok(NULL, ",");
				i++;
			}
		}
		struct idx_head head;
		TAILQ_INIT(&head);
		snp_data snp_info = read_snp_file(snp_file);
		get_multiple_snp_idx(&snp_info, snp_names, n_elems, &head, NULL);
		
		struct idx_node* in;
		if(strcmp(idx_str, "NA") == 0) {
			TEST_ASSERT_TRUE(TAILQ_EMPTY(&head));
		} else {
			i = 0;
			TAILQ_FOREACH(in, &head, nodes) {
				TEST_ASSERT_EQUAL_UINT(exp_idx[i], in->idx);
				i++;
			}
		}
		free_snp_data(&snp_info);
		free_idx_list(&head);
	}
	free(buf);
	fclose(fp);
}

void test_multiple_index_ind(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/multiple_index_test/ind", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char* ind_file;
		size_t n_elems;
		char* ind_str;
		char* pop_str;
		char* idx_str;
		size_t cur_beg = 0;
		int n_col = 0;
		// split string
		for(size_t c_i = 0; c_i < bufsize; c_i++) {
			if((buf[c_i] == '\t') || buf[c_i] == '\n') {
				buf[c_i] = '\0';
				switch(n_col) {
					case 0:
						ind_file = &buf[cur_beg];
						break;
					case 1:
						n_elems = atoi(&buf[cur_beg]);
						break;
					case 2:
						ind_str = &buf[cur_beg];
						break;
					case 3:
						pop_str = &buf[cur_beg];
						break;
					case 4:
						 idx_str = &buf[cur_beg];
						 break;
				}
				cur_beg = c_i + 1;
				n_col++;
			}
		}

		char** ind_names = malloc(sizeof(char*) * n_elems);
		char** pop_names = malloc(sizeof(char*) * n_elems);
		size_t* exp_idx = malloc(sizeof(size_t) * n_elems);
	
		char* elem = strtok(ind_str, ",");
		size_t i = 0;
		while(elem) {
			ind_names[i] = strdup(elem);
			elem = strtok(NULL, ",");
			i++;
		}

		elem = strtok(pop_str, ",");
		i = 0;
		while(elem) {
			pop_names[i] = strdup(elem);
			elem = strtok(NULL, ",");
			i++;
		}

		if(strcmp(idx_str, "NA") != 0) {
			elem = strtok(idx_str, ",");
			i = 0;
			while(elem) {
				exp_idx[i] = atoi(elem);
				elem = strtok(NULL, ",");
				i++;
			}
		}
		
		struct idx_head head;
		TAILQ_INIT(&head);
		ind_data ind_info = read_ind_file(ind_file);
		get_multiple_ind_idx(&ind_info, ind_names, pop_names, n_elems, &head, NULL);
		struct idx_node* in;
		if(strcmp(idx_str, "NA") == 0) {
			TEST_ASSERT_TRUE(TAILQ_EMPTY(&head));
		} else {
			i = 0;
			TAILQ_FOREACH(in, &head, nodes) {
				TEST_ASSERT_EQUAL_UINT(exp_idx[i], in->idx);
				i++;
			}
		}
		free_ind_data(&ind_info);
		free_idx_list(&head);
	}
	free(buf);
	fclose(fp);
}

void test_multiple_index_in_ind_file(void) {
	for(int i = 1; i <= 10; i++) {
		char ind_file[46];
		sprintf(ind_file, "test/data/multidim_ind/multi_index_ind_%d.ind", i);
		ind_data ind_info = read_ind_file(ind_file);
		TEST_ASSERT_TRUE(ind_info.length > 0);
		free_ind_data(&ind_info);
	}
}

void test_append_ind(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/append_ind/append_table.tsv", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char* file1 = strtok(buf, "\t\n");
		char* file2 = strtok(NULL, "\t\n");
		char* file_exp = strtok(NULL, "\t\n");

		ind_data ind1 = read_ind_file(file1);
		ind_data ind2 = read_ind_file(file2);
		ind_data ind_exp = read_ind_file(file_exp);
		ind_data ind_actual;

		append_ind_data(&ind1, &ind2, &ind_actual);
		TEST_ASSERT_EQUAL_UINT(ind_exp.length, ind_actual.length);
		TEST_ASSERT_EQUAL_UINT32(ind_exp.hash, ind_actual.hash);
		TEST_ASSERT_EQUAL_STRING_ARRAY(ind_exp.ind_id, ind_actual.ind_id, ind_actual.length);
		TEST_ASSERT_EQUAL_STRING_ARRAY(ind_exp.sex, ind_actual.sex, ind_actual.length);
		TEST_ASSERT_EQUAL_STRING_ARRAY(ind_exp.population, ind_actual.population, ind_actual.length);

		free_ind_data(&ind1);
		free_ind_data(&ind2);
		free_ind_data(&ind_exp);
		free_ind_data(&ind_actual);
	}
}

void test_intersect_snp(void) {
	char* buf = NULL;
	size_t bufsize = 0;
	FILE* fp = fopen("test/data/snp_isct/snp_isct.tsv", "r");
	while(getline(&buf, &bufsize, fp) != -1) {
		char* file1;
		char* file2;
		char* idx_str1;
		char* idx_str2;
		size_t cur_beg = 0;
		int n_col = 0;
		for(size_t c_i = 0; c_i < bufsize; c_i++) {
			if((buf[c_i] == '\t') || buf[c_i] == '\n') {
				buf[c_i] = '\0';
				switch(n_col) {
					case 0:
						file1 = &buf[cur_beg];
						break;
					case 1:
						file2 = &buf[cur_beg];
						break;
					case 2:
						idx_str1 = &buf[cur_beg];
						break;
					case 3:
						idx_str2 = &buf[cur_beg];
						break;
				}
				cur_beg = c_i + 1;
				n_col++;
			}
		}

		snp_data snp1 = read_snp_file(file1);
		snp_data snp2 = read_snp_file(file2);
		struct idx_head head_idx1;
		struct idx_head head_idx2;
		TAILQ_INIT(&head_idx1);
		TAILQ_INIT(&head_idx2);
		intersect_snp_data(&snp1, &snp2, &head_idx1, &head_idx2);
		struct idx_node* in;
		size_t i = 0;
		TAILQ_FOREACH(in, &head_idx1, nodes) {
			size_t elem_exp;
			if(i == 0) {
				elem_exp = atoi(strtok(idx_str1, ","));
			} else {
				elem_exp = atoi(strtok(NULL, ","));
			}
			TEST_ASSERT_EQUAL_UINT(elem_exp, in->idx);
			i++;
		}

		free_snp_data(&snp1);
		free_snp_data(&snp2);
		free_idx_list(&head_idx1);
		free_idx_list(&head_idx2);
	}
}

int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_equality_egn_pam);
	RUN_TEST(test_goto_equality_pam_egn);
	RUN_TEST(test_ind_hash);
	RUN_TEST(test_snp_hash);
	RUN_TEST(test_ind_filtered_hash);
	RUN_TEST(test_snp_filtered_hash);
	RUN_TEST(test_read_write_eqv_pam);
	RUN_TEST(test_read_write_eqv_tgn);
	RUN_TEST(test_read_write_eqv_egn);
	RUN_TEST(test_intersect_snp);
	RUN_TEST(test_get_idx_snp);
	RUN_TEST(test_get_idx_ind);
	RUN_TEST(test_multiple_index_snp);
	RUN_TEST(test_multiple_index_ind);
	RUN_TEST(test_multiple_index_in_ind_file);
	RUN_TEST(test_append_ind);
	RUN_TEST(test_pop_filtered_hash);
	RUN_TEST(test_chr_filtered_hash);
	RUN_TEST(test_range_filtered_hash);
    return UNITY_END();
}
