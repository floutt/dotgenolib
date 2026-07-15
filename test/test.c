#include "unity.h"
#include "dotgeno.h"

void setUp(void) {
    // set stuff up here
}

void tearDown(void) {
    // clean stuff up here
}

void test_equality_egn_pam(void) {
	char* file_base_egn = "../test/data/PAM_EGN/egn_example_%i.%s";
	char* file_base_pam = "../test/data/PAM_EGN/pam_example_%i.%s";
	char buf[41];
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
					TEST_FAIL();
				}
				break;
			}
			TEST_ASSERT_EQUAL_UINT8_ARRAY(pam_record, egn_record, ind_info_pam.length);
		}
		free_snp_data(&snp_info_pam);
		free_snp_data(&snp_info_egn);
		free_ind_data(&ind_info_pam);
		free_ind_data(&ind_info_egn);
		close_pam_file_reader(&pfr);
		close_egn_file_reader(&efr);
	}
}

// not needed when using generate_test_runner.rb
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_equality_egn_pam);
    return UNITY_END();
}
