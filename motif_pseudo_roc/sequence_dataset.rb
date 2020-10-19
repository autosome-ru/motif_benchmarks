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
end
