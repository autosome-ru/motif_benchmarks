class FastaSequence
  attr_reader :sequence, :header
  def initialize(sequence, header)
    @sequence = sequence.upcase
    @header = header
  end

  def length
    sequence.length
  end

  # yields nucleotide and weight for each position
  def each_position(&block)
    sequence.each_char(&block)
  end

  def self.each_in_file(filename)
    return enum_for(:each_in_file, filename)  unless block_given?
    File.open(filename) do |f|
      f.each_line.lazy.map(&:strip).slice_before{|line|
        line.start_with?('>')
      }.each do |lines|
        header, *sequence_parts = lines
        yield self.new(sequence_parts.join, header)
      end
    end
  end
end
