#ifndef DOTGENO_H
#define DOTGENO_H

#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <sys/queue.h>

#ifdef __cplusplus
extern "C" {
#endif

/*! @brief NAN value encoding */
#define NAN_VAL 3

/*! @typedef
 * @brief Struct for indexing ind_data objects
 *
 * This is used to access information from ".ind" files.
 */
typedef struct {
	char* ind_id; /**< Intrapopulation identifier for individual */
	char* ind_pop; /**< Unique population identifier*/
} ind_idx;


struct idx_node {
	size_t idx;
	STAILQ_ENTRY(idx_node) nodes;
};

/*
 * @brief queue.h linked list head for ind_data and snp_data indices
 */
STAILQ_HEAD(idx_head, idx_node);

struct str_node {
	char* str;
	STAILQ_ENTRY(str_node) nodes;
};

/*
 * @brief queue.h linked list head for strings
 */
STAILQ_HEAD(str_list_head, str_node);

struct ind_idx_node {
	ind_idx* iidx;
	STAILQ_ENTRY(ind_idx_node) nodes;
};

/*
 * @brief queue.h linked list head for storing ind_idx objects
 */
STAILQ_HEAD(ind_idx_head, ind_idx_node);

/*! @brief string to size_t hash table for use with snp_data structs */
typedef struct id_map_str id_map_str;

/*! @brief ind_idx to size_t hash table for use with ind_data structs */ 
typedef struct id_map_ind id_map_ind;

/*! @typedef
 * @brief struct containing information encoded in ".snp" files
 */
typedef struct {
	size_t length; /**< Total number of elements for each array in struct */
	char** var_id; /**< Array of variant identifiers */
	char** chr; /**< Array of chromosomes of variants */
	double* pos_morgans; /**< Array of genetic positions in Morgans */
	uint64_t* pos; /**< Array of genetic base pair positions */
	char** ref; /**< Array of reference alleles */
	char** alt; /**< Array of alternative alleles */ 
	id_map_str* rev_idx; /**< Hash table which returns index of the given var_id */
	uint32_t hash; /**< Hash value for struct */
} snp_data;

/*! @typedef
 * @brief struct containing information encoded in ".ind" files
 */
typedef struct {
	size_t length; /**< Total number of elements for each array in struct */
	char** ind_id; /**< Array of intrapopulation identifiers for individuals */
	char** sex; /**< Array of sex for individuals */
	char** population; /**< Array of sex for individuals */
	id_map_ind* rev_idx; /**< Hash table which returns index of the given ind_idx struct */
	uint32_t hash; /**< Hash value for struct */
} ind_data;

/*! @typedef
 * @brief struct containing information necessary to read PACKEDANCESTRYMAP files
 */
typedef struct {
	bool is_hdr_read; /**< Whether the header of the PACKEDANCESTRYMAP file has been read */
	bool is_open; /**< Whether the PACKEDANCESTRYMAP file has been opened */
	size_t idx; /**< Index of the current record being iterated */
	size_t record_size; /**< Size of each PACKEDANCESTRYMAP record */
	size_t n_ind; /**< Number of individuals in the file */
	size_t n_snp; /**< Number of genetic variants in the file */
	FILE* fp;  /**< File pointer to PACKEDANCESTRYMAP file */ 
} pam_file_reader;

/*! @typedef
 * @brief struct containing information necessary to read EIGENSTRAT files
 */
typedef struct {
	bool is_open; /**< Whether the EIGENSTRAT file has been opened */
	size_t idx; /**< Index of the current record being iterated */
	size_t n_ind; /**< Number of individuals in the file */
	size_t n_snp; /**< Number of genetic variants in the file */  
	FILE* fp; /**< File pointer to EIGENSTRAT file */
} egn_file_reader;

typedef struct {
	bool hdr_written; /**< Whether the header of the PACKEDANCESTRYMAP file has been written */
	bool is_open; /**< Whether the PACKEDANCESTRYMAP file has been opened */
	size_t n_ind; /**< Number of individuals to write in the file */ 
	size_t n_snp; /**< Number of genetic variants to write in the file */ 
	size_t n_written_snp; /**< Number of variants written to the file */
	size_t record_size; /**< Size of each PACKEDANCESTRYMAP record */
	FILE* fp; /**< File pointer to PACKEDANCESTRYMAP file */
} pam_file_writer;

typedef struct {
	bool is_open; /**< Whether the EIGENSTRAT file has been opened */
	size_t n_ind; /**< Number of individuals to write in the file */ 
	size_t n_snp; /**< Number of genetic variants to write in the file */
	size_t n_written_snp; /**< Number of variants written to the file */
	FILE* fp; /**< File pointer to EIGENSTRAT file */
} egn_file_writer;

/**
 * @brief writes the header record of PACKEDANCESTRYMAP file
 *
 * @param[out] pfw pam_file_writer object to write header to
 * @param[in] snp_info corresponding variant information for pfw
 * @param[in] ind_info corresponding individual information for pfw
 */
void write_pam_header(pam_file_writer* pfw, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief writes variant record to PACKEDANCESTRYMAP file 
 *
 * @param[out] pfw pam_file_writer object to write record to
 * @param[in] dosages array of dosages. Should be of length pfw->n_ind
 */
void write_pam_record(pam_file_writer* pfw, uint8_t* dosages);

/**
 * @brief writes variant record to EIGENSTRAT file 
 *
 * @param[out] efw egn_file_writer object to write record to
 * @param[in] dosages array of dosages. Should be of length efw->n_ind
 */
void write_egn_record(egn_file_writer* efw, uint8_t* dosages);

/*! @typedef
 * @brief struct containing information encoded in PACKEDANCESTRYMAP file headers
 */
typedef struct {
	size_t n_ind; /**< Number of individuals in the file */
	size_t n_snp; /**< Number of genetic variants in the file */ 
	uint32_t ind_hash; /**< Expected hash value for corresponding ".ind" file */
	uint32_t snp_hash; /**< Expected hash value for corresponding ".snp" file */
} hdr_data;

/**
 *  @brief reads information from .snp file and stores it in snp_data struct
 *
 *  @param[in] filename path to .snp file
 *
 *  @return snp_data struct containing information from .snp file
 */ 
snp_data read_snp_file(char* filename);

/**
 *  @brief reads information from .ind file and stores it in ind_data struct
 *
 *  @param[in] filename path to .ind file
 *
 *  @return ind_data struct containing information from .ind file
 */ 
ind_data read_ind_file(char* filename);

/**
 * @brief writes contents of snp_data object to file
 *
 * @param[in] snp_info pointer to snp_data object
 * @param[in] filename path to output file
 *
 * @return status code indicating whether or not file was successfully read and closed
 * @retval 0 not successfully closed 
 * @retval EOF file not successfully closed
 */
short write_snp_data(snp_data* snp_info, char* filename);

/**
 * @brief writes contents of ind_data object to file
 *
 * @param[in] ind_info pointer to ind_data object
 * @param[in] filename path to output file
 *
 * @return status code indicating whether or not file was successfully read and closed
 * @retval 0 not successfully closed 
 * @retval EOF file not successfully closed
 */
short write_ind_data(ind_data* ind_info, char* filename);

/**
 * @brief gets the index for snp with the identifier var_name
 *
 * @param[in] snp_info pointer to snp_data struct storing variant information
 * @param[in] var_name ID of desired variant
 * @param[out] idx pointer to where index will be stored
 *
 * @return status code indicating whether or not variant was found in snp_info
 * @retval -1 variant not found
 * @retval 0 variant found
 */
short get_snp_idx(snp_data* snp_info, char* var_name, size_t* idx);

/**
 * @brief gets the index for an individual with with certain individual ID AND population
 *
 * @param[in] ind_info pointer to ind_data struct storing variant information
 * @param[in] ind_id intrapopulation identifier for individual
 * @param[in] ind_pop population of individual
 * @param[out] idx pointer to where index will be stored
 *
 * @return status code indicating whether or not individual was found in ind_info
 * @retval -1 individual not found
 * @retval 0 individual found
 */
short get_ind_idx(ind_data* ind_info, char* ind_id, char* ind_pop, size_t* idx);

/**
 * @brief returns the numerical indices of a set of var_names
 * 
 * This function will only return indices for elements of var_names found in the snp_data object
 *
 * @param[in] snp_info snp_data object to be queried
 * @param[in] var_names array of strings representing variable names
 * @param[in] length length of the var_names array
 * @param[out] head_idx head of the index linked list where indexes will be stored
 * @param[out] head_str head of string linked list where variant names not found will be stored. Set to NULL if you don't want to store this information
 */
void get_multiple_snp_idx(snp_data* snp_info, char** var_names, size_t length, struct idx_head* head_idx, struct str_list_head* head_str);

/**
 * @brief returns the numerical indices of a set of individuals
 * 
 * This function will only return indices for elements of individuals found in the ind_data object
 *
 * @param[in] ind_info ind_data object to be queried
 * @param[in] ind_ids array of strings representing the individual IDs
 * @param[in] ind_pops array of strings representing the individuals' respective populations 
 * @param[in] length length of the ind_ids and ind_pops arrays
 * @param[out] head_idx head of the index linked list where indexes will be stored
 * @param[out] head_iidx head of ind_idx linked list where individual identifier structs of individuals not found will be stored. Set to NULL if you don't want to store this information
 */
void get_multiple_ind_idx(ind_data* ind_info, char** ind_ids, char** ind_pops, size_t length, struct idx_head* head_idx, struct ind_idx_head* head_iidx);

/**
 * @brief filters snp_data object to only those elements found in the provided list
 *
 * @param[in] snp_in snp_data object to be filtered
 * @param[out] snp_out snp_data object in which to write the filtered data to
 * @param[in] head head of linked list containing the indices from snp_in to include in snp_out
 *
 * @return status code indicating whether snp_out was modified
 * @retval 0 snp_out successfully modified
 * @retval -1 snp_out not modified, linked list was of length 0
 */
short filter_snp_data(snp_data* snp_in, snp_data* snp_out, struct idx_head* head);

/**
 * @brief filters ind_data object to only those elements found in the provided list
 *
 * @param[in] ind_in ind_data object to be filtered
 * @param[out] ind_out ind_data object in which to write the filtered data to
 * @param[in] head head of linked list containing the indices from ind_in to include in ind_out
 *
 * @return status code indicating whether ind_out was modified
 * @retval 0 ind_out successfully modified
 * @retval -1 ind_out not modified, linked list was of length 0
 */
short filter_ind_data(ind_data* ind_in, ind_data* ind_out, struct idx_head* head);

/**
 * @brief finds the intersecting rows of two snp_data objects
 *
 * @param[in] snp1 snp_data object 1
 * @param[in] snp2 snp_data object 2
 * @param[out] head1 head of the linked list where snp1 indices of intersecting rows will be stored
 * @param[out] head2 head of the linked list where snp2 indices of intersecting rows will be stored
 *
 * @return number of intersecting rows between snp1 and snp2
 */
size_t intersect_snp_data(snp_data* snp1, snp_data* snp2, struct idx_head* head1, struct idx_head* head2);

/**
 * @brief appends two ind_data objects
 *
 * @param[in] ind1 first ind_data object
 * @param[in] ind2 second ind_data object
 * @param[out] ind_out ind_data object whete ind1 + ind2 will be stored
 */
void append_ind_data(ind_data* ind1, ind_data* ind2, ind_data* ind_out);

/**
 * @brief opens PACKEDANCESTRYMAP file and initializes reader
 *
 * @param[in] filename path to PACKEDANCESTRYMAP file
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] ind_info pointer to corresponding ind_data object
 *
 * @return object with file pointer and information necessary to read PACKEDANCESTRYMAP file 
 */
pam_file_reader pam_file_reader_init(char* filename, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief opens EIGENSTRAT file and initializes reader
 *
 * @param[in] filename path to EIGENSTRAT file
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] ind_info pointer to corresponding ind_data object
 *
 * @return object with file pointer and information necessary to read EIGENSTRAT file 
 */
egn_file_reader egn_file_reader_init(char* filename, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief initializes PACKEDANCESTRYMAP file writer
 *
 * @param[in] filename path to PACKEDANCESTRYMAP file
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] ind_info pointer to corresponding ind_data object
 *
 * @return object with file pointer and information necessary to read PACKEDANCESTRYMAP file 
 */
pam_file_writer pam_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief initializes EIGENSTRAT file writer
 *
 * @param[in] filename path to EIGENSTRAT file
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] ind_info pointer to corresponding ind_data object
 *
 * @return object with file pointer and information necessary to read EIGENSTRAT file 
 */
egn_file_writer egn_file_writer_init(char* filename, snp_data* snp_info, ind_data* ind_info);

/**
 * @brief closes pam_file_reader object
 *
 * @param[in,out] pf pointer to pam_file_reader object
 *
 * @return status code indicating whether or not the file was successfully closed
 * @retval 0 file successfully closed
 * @retval EOF file not successfully closed
 */
int close_pam_file_reader(pam_file_reader* pf);

/**
 * @brief closes egn_file_reader object
 *
 * @param[in,out] ef pointer to egn_file_reader object
 *
 * @return status code indicating whether or not the file was successfully closed
 * @retval 0 file successfully closed
 * @retval EOF file not successfully closed
 */
int close_egn_file_reader(egn_file_reader* ef);

/**
 * @brief closes pam_file_writer object
 *
 * @param[in,out] pfw pointer to pam_file_writer object
 *
 * @return status code indicating whether or not the file was successfully closed
 * @retval 0 file successfully closed
 * @retval EOF file not successfully closed
 */
int close_pam_file_writer(pam_file_writer* pfw);

/**
 * @brief closes egn_file_writer object
 *
 * @param[in,out] efw pointer to egn_file_writer object
 *
 * @return status code indicating whether or not the file was successfully closed
 * @retval 0 file successfully closed
 * @retval EOF file not successfully closed
 */
int close_egn_file_writer(egn_file_writer* efw);

/**
 * @brief reads header record of PACKEDANCESTRY map file.
 * This must be read before any other record is read
 *
 * @param[in,out] pf pointer to pam_file_reader object 
 *
 * @return hdr_data object containing important header information
 */
hdr_data read_pam_header(pam_file_reader* pf);

/**
 * @brief reads record of dosages of the current genetic variant from pam_file_reader object and iterates to the next one
 *
 * @param[in,out] pf pointer to pam_file_reader object
 *
 * @return array of allelic dosages of length pf->length. Dosages must be either 0,1,2, or NAN_VAL.
 */
uint8_t* read_pam_record(pam_file_reader* pf);

/**
 * @brief reads record of dosages of the current genetic variant from egn_file_reader object and iterates to the next one
 *
 * @param[in,out] ef pointer to egn_file_reader object
 *
 * @return array of allelic dosages of length ef->length. Dosages must be either 0,1,2, or NAN_VAL.
 */
uint8_t* read_egn_record(egn_file_reader* ef);

/**
 * @brief sends egn_file_reader object to the location of a particular variant
 *
 * @param[in,out] ef pointer to egn_file_reader object
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] var_name variant to transport object to
 * @return status code indicating whether or not egn_file_reader was relocated to variant var_name
 * @retval -1 variant not found
 * @retval 0 successful relocation to variant var_name
 */
short goto_var_egn(egn_file_reader* ef, snp_data* snp_info, char* var_name);

/**
 * @brief sends pam_file_reader object to the location of a particular variant
 *
 * @param[in,out] pf pointer to pam_file_reader object
 * @param[in] snp_info pointer to corresponding snp_data object
 * @param[in] var_name variant to transport object to
 * @return status code indicating whether or not pam_file_reader was relocated to variant var_name
 * @retval -1 variant not found
 * @retval 0 successful relocation to variant var_name
 */
short goto_var_pam(pam_file_reader* pf, snp_data* snp_info, char* var_name);

/**
 * @brief frees snp_data object
 */
void free_snp_data(snp_data* snp_info);

/**
 * @brief frees ind_data object
 */
void free_ind_data(ind_data* ind_info);

/**
 * @brief empties and frees index linked list
 */
void free_idx_list(struct idx_head* head);

/**
 * @brief emptiies and frees string linked list
 */
void free_str_list(struct str_list_head* head);

/**
 * @brief empties and frees ind_idx linked list
 */
void free_ind_idx_list(struct ind_idx_head* head);

#ifdef __cplusplus
}
#endif

#endif
