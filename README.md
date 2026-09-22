# dotgenolib

`dotgenolib` is a fast, minimalist C library for reading and writing PACKEDANCESTRYMAP and EIGENSTRAT files.

## Installation
Installation of this software is straightforward. No external dependencies dependencies are required outside of the C standard library and POSIX standard. Simply run the following commands:
```sh
make
make test  # run unit tests
make install
make clean
```

## Overview of library API
### Read functions
| Operation | EIGENSTRAT | PACKEDANCESTRYMAP | TGENO |
|---|---|---|---|
| Initialize reader | `egn_file_reader_init` | `pam_file_reader_init` | `tgn_file_reader_init` |
| Read header | N/A | `read_pam_header` | `read_tgn_header` |
| Read record | `read_egn_record` | `read_pam_record` | `read_tgn_record` |
| Free record | `free` | `free` | `free` |
| Random access | `goto_var_egn` | `goto_var_pam` | `goto_ind_tgn` |
| Close reader | `close_egn_file_reader` | `close_pam_file_reader` | `close_tgn_file_reader` |

### Write functions
| Operation | EIGENSTRAT | PACKEDANCESTRYMAP | TGENO |
|---|---|---|---|
| Initialize writer | `egn_file_writer_init` | `pam_file_writer_init` | `tgn_file_writer_init` |
| Write header | N/A | `write_pam_header` | `write_tgn_header` |
| Write record | `write_egn_record` | `write_pam_record` | `write_tgn_record` |
| Close writer | `close_egn_file_writer` | `close_pam_file_writer` | `close_tgn_file_writer` |

### Metadata files

| Property | `.snp` | `.ind` |
|---|---|---|
| Contains | Genetic variant metadata | Individual metadata |
| Rows represent | Genetic variants | Individuals |
| Read function | `read_snp_file` | `read_ind_file` |
| Write function | `write_snp_data` | `write_ind_data` |
| Free function | `free_snp_data` | `free_ind_data` |

### Record structure

| Property | EIGENSTRAT | PACKEDANCESTRYMAP | TGENO |
|---|---|---|---|
| Record represents | Genetic variants | Genetic variants | Individuals |
| Records | `n_snp` | `n_snp` | `n_ind` |
| Values per record | `n_ind` | `n_ind` | `n_snp` |
| Storage | Text | 2-bit packed | 2-bit packed |

### Filtering / Index functions

#### Lookups (single item)

| Function | Description |
| --- | --- |
| `get_snp_idx` | Looks up the index of a single variant by its ID. |
| `get_ind_idx` | Looks up the index of a single individual by its individual and population IDs. |

#### Lookups (multiple items)

| Function | Description |
| --- | --- |
| `get_multiple_snp_idx` | Returns the indices for a list of variant IDs. Matches are appended to `head_idx`; any names not found are optionally collected in `head_str`. |
| `get_multiple_ind_idx` | Returns the indices for a list of (ID, population) pairs. Matches are appended to `head_idx`; unmatched individuals are optionally collected in `head_iidx`. |
| `get_multiple_pops` | Returns the indices of every individual belonging to any population in `ind_pops`. Populations with no matching individuals are optionally collected in `head_nopop`. |
| `get_multiple_sex` | Returns the indices of every individual with the given sex. |
| `get_multiple_chrs` | Returns the indices of every variant located on any chromosome in `chrs`. |
| `get_multiple_ranges` | Returns the indices of every variant falling within any of the given (chromosome, start, end) ranges. Ranges are inclusive on both ends. |

#### Applying filters

| Function | Description |
| --- | --- |
| `filter_snp_data` | Writes to `snp_out` only the rows of `snp_in` whose indices appear in `head`. Returns `-1` (unmodified) if `head` is empty, `0` on success. |
| `filter_ind_data` | Writes to `ind_out` only the rows of `ind_in` whose indices appear in `head`. Returns `-1` (unmodified) if `head` is empty, `0` on success. |
| `intersect_snp_data` | Finds variants common to `snp1` and `snp2`, storing the matching indices for each in `head1`/`head2` respectively. Returns the number of intersecting rows. |

#### Cleanup

| Function | Frees |
| --- | --- |
| `free_idx_list` | An `idx_head` linked list of numeric indices (the result of any `get_*`/`get_multiple_*` lookup). |
| `free_str_list` | A `str_list_head` linked list of strings (e.g. unmatched variant IDs or populations). |
| `free_ind_idx_list` | An `ind_idx_head` linked list of `ind_idx` structs (e.g. unmatched individuals). |
