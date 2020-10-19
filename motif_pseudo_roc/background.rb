require_relative 'support'

def local_mono_background(sequence_dataset)
  counts = Hash.new(0)
  sequence_dataset.each_sequence{|fasta_sequence|
    fasta_sequence.each_position{|letter|
      counts[letter] += 1
    }
  }
  counts = counts.select{|kmer, cnt| Nucleotides.include?(kmer) }
  normed_hash(counts_symmetrized(counts)).values_at(*Nucleotides)
end

def local_di_background(sequence_dataset)
  counts = Hash.new(0)
  sequence_dataset.each_sequence{|fasta_sequence|
    fasta_sequence.each_position.each_cons(2){|letter_1, letter_2|
      diletter = "#{letter_1}#{letter_2}"
      counts[diletter] += 1
    }
  }
  counts = counts.select{|kmer, cnt| Dinucleotides.include?(kmer) }
  normed_hash(counts_symmetrized(counts)).values_at(*Dinucleotides)
end
