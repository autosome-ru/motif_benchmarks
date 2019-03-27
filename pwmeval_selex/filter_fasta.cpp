#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cctype>
#include <cstring>

bool has_only_acgt(const std::string& seq) {
	for (size_t i = 0; i < seq.length(); ++i){
		char letter = toupper(seq[i]);
		if (!(letter == 'A' || letter == 'C' || letter == 'G' || letter == 'T')) {
			return false;
		}
	}
	return true;
}

void filter_fasta(std::istream& input, std::ostream& output, bool only_acgt = false, size_t seq_length = 0) {
	bool skip = false;
	std::string seq_id, seq;
	while (input.good()) {
		std::string line;
		std::getline(input, line);
		if (line.length() > 0) {
			if (line[0] == '>') {
				if (seq_id.length() > 0 || seq.length() > 0) {
					skip = (only_acgt && !has_only_acgt(seq)) || ((seq_length != 0) && (seq.length() != seq_length));
					if (!skip) {
						output << seq_id << std::endl << seq << std::endl;
					}
					seq.clear();
				}
				seq_id = line;
			} else {
				seq += line;
			}
		} 
	}
  if (seq_id.length() > 0 || seq.length() > 0) {
    skip = (only_acgt && !has_only_acgt(seq)) || ((seq_length != 0) && (seq.length() != seq_length));
    if (!skip) {
      output << seq_id << std::endl << seq << std::endl;
    }
  }
}

int main(int argc, char **argv) {
  size_t seq_length;
  bool only_acgt;
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <filename or - for stdin> <sequence length = integer|no> <only acgt = yes|no>" << std::endl;
    exit(1);
  }
  
  if (!strcmp(argv[2], "no")) {
    seq_length = 0;
  } else {
    seq_length = atoi(argv[2]);
  }

  if (!strcmp(argv[3], "yes")) {
    only_acgt = true;
  } else if (!strcmp(argv[3], "no")) {
    only_acgt = false;
  } else {
    std::cerr << "Usage: " << argv[0] << " <filename or - for stdin> <sequence length = integer|no> <only acgt = yes|no>" << std::endl;
    exit(1); 
  }
	if (!strcmp(argv[1], "-")) {
		filter_fasta(std::cin, std::cout, only_acgt, seq_length);
	} else {
		std::ifstream fasta_file(argv[1], std::ifstream::in);
		if (fasta_file.fail()) {
			std::cerr << "Failed to open file" << std::endl;
			exit(1);
		}
		filter_fasta(fasta_file, std::cout, only_acgt, seq_length);
	}
  return 0;
}
