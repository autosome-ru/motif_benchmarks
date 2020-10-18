require_relative 'support'
require_relative 'fasta_sequence'

class SequenceDataset
  attr_reader :filename

  def initialize(filename)
    @filename = filename
  end

  def each_sequence
    return enum_for(:each_sequence)  unless block_given?
    FastaSequence.each_in_file(filename){|fasta_sequence|
      yield fasta_sequence
    }
  end

  def local_mono_background
    @local_mono_background ||= begin
      counts = Hash.new(0)
      each_sequence{|fasta_sequence|
        fasta_sequence.each_position{|letter|
          counts[letter] += 1
        }
      }
      normed_hash(counts_symmetrized(counts)).values_at(*Nucleotides)
    end
  end

  def local_di_background
    @local_di_background ||= begin
      counts = Hash.new(0)
      each_sequence{|fasta_sequence|
        fasta_sequence.each_position.each_cons(2){|letter_1, letter_2|
          diletter = "#{letter_1}#{letter_2}"
          counts[diletter] += 1
        }
      }
      normed_hash(counts_symmetrized(counts)).values_at(*Dinucleotides)
    end
  end
end
