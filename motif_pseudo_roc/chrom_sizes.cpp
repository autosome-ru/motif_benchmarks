#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cctype>
#include <cstring>

class ContigInfo {
public:
  std::string contig;
  int length;
  ContigInfo(std::string contig, int length) : contig(contig), length(length) {  }
  friend std::ostream& operator<<(std::ostream& out, const ContigInfo& info);
};

std::ostream& operator<<(std::ostream& out, const ContigInfo& info) {
  out << info.contig << "\t" << info.length << std::endl;
  return out;
}

std::vector<ContigInfo> count_fasta_sizes(std::istream& input) {
  std::vector<ContigInfo> result;
  std::string seq_id;
  int seq_len = 0;
  while (input.good()) {
    std::string line;
    std::getline(input, line);
    if (line.length() > 0) {
      if (line[0] == '>') {
        if (seq_id.length() > 0) {
          result.push_back(ContigInfo(seq_id, seq_len));
        }
        seq_id = line.substr(1);
        seq_len = 0;
      } else {
        seq_len += line.length();
      }
    }
  }
  if (seq_id.length() > 0) {
    result.push_back(ContigInfo(seq_id, seq_len));
  }
  return result;
}

void print_sizes(const std::vector<ContigInfo>& contig_sizes, std::ostream& output) {
  for (auto iter = contig_sizes.begin(); iter != contig_sizes.end(); ++iter) {
    output << *iter;
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <filename or - for stdin>" << std::endl;
    exit(1);
  }

  if (!strcmp(argv[1], "-")) {
    print_sizes(count_fasta_sizes(std::cin), std::cout);
  } else {
    std::ifstream fasta_file(argv[1], std::ifstream::in);
    if (fasta_file.fail()) {
      std::cerr << "Failed to open file" << std::endl;
      exit(1);
    }
    print_sizes(count_fasta_sizes(fasta_file), std::cout);
  }
  return 0;
}
