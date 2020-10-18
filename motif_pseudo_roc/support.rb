Nucleotides = %w[A C G T]
Dinucleotides = Nucleotides.product(Nucleotides).map(&:join)

def counts_symmetrized(counts)
  result = counts.dup.tap{|h| h.default = 0 }
  counts.each{|kmer, cnt|
    revcomp_kmer = kmer.tr('ACGT','TGCA').reverse
    result[revcomp_kmer] += cnt
  }
  result
end

def normed_hash(counts)
  total_count = counts.map{|kmer, cnt| cnt }.sum(0.0)
  counts.map{|kmer, cnt| [kmer, cnt.to_f / total_count] }.to_h
end
