#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cctype>
#include <cstring>

void count_fasta_sizes(std::istream& input, std::ostream& output) {
  std::string seq_id;
  int seq_len = 0;
  while (input.good()) {
    std::string line;
    std::getline(input, line);
    if (line.length() > 0) {
      if (line[0] == '>') {
        if (seq_id.length() > 0) {
          output << seq_id << "\t" << seq_len << std::endl;
        }
        seq_id = line.substr(1);
        seq_len = 0;
      } else {
        seq_len += line.length();
      }
    }
  }
  if (seq_id.length() > 0) {
    output << seq_id << "\t" << seq_len << std::endl;
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <filename or - for stdin>" << std::endl;
    exit(1);
  }

  if (!strcmp(argv[1], "-")) {
    count_fasta_sizes(std::cin, std::cout);
  } else {
    std::ifstream fasta_file(argv[1], std::ifstream::in);
    if (fasta_file.fail()) {
      std::cerr << "Failed to open file" << std::endl;
      exit(1);
    }
    count_fasta_sizes(fasta_file, std::cout);
  }
  return 0;
}
