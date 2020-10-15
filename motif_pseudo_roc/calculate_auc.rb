####  Input:
####    >seqHeader:seqLength
####    pvalue <TAB> any other info (e.g. position and orientation)
####    >seqHeader:seqLength
####    pvalue	position	orientation
####  `seqHeader` has no constraints (such as uniqueness)
#### Example:
## >ANDR_HUMAN.PEAKS033090.pics.control:100
## 0.001       5       -
## >ANDR_HUMAN.PEAKS033090.pics.control:113
## 3.5e-4       173     +

require 'json'

def median(arr)
  raise 'Median of an empty array is undefined'  if arr.empty?
  sorted_arr = arr.sort
  if sorted_arr.size.odd?
    center_index = sorted_arr.size / 2
    sorted_arr[center_index]
  else
    center_index_ceiled = sorted_arr.size / 2
    (sorted_arr[center_index_ceiled] + sorted_arr[center_index_ceiled - 1]) / 2.0
  end
end

def correct_pvalues(pvalues, median_sequence_length:, model_length:)
  pvalues.map{|pvalue|
    1.0 - (1.0 - pvalue) ** (2 * (median_sequence_length - model_length + 1))
  }
end

def calculate_auc(points, xy_keys: [0, 1])
  xy_coords = points.map{|point| point.values_at(*xy_keys) }
  xy_coords.sort_by{|x,y|
    x
  }.each_cons(2).map{|(x1, y1), (x2, y2)|
    (x2 - x1) * (y1 + y2) / 2.0
  }.inject(0.0, &:+)
end

def roc_curve(pvalues)
  num_positive_sequences = pvalues.size # all sequences in dataset are positive

  roc_points = []
  roc_points << {tpr: 0.0, fpr: 0.0}
  # list pvalues from the most stringent threshold (no one sequence taken as a positive one) to the weakest
  sorted_pvalues = pvalues.sort
  sorted_pvalues.each_with_index{|corrected_pvalue, index|
    num_positive_sequences_taken = index + 1
    tpr = num_positive_sequences_taken.to_f / num_positive_sequences
    fpr = corrected_pvalue
    roc_points << {tpr: tpr, fpr: fpr}
  }
  roc_points << {tpr: 1.0, fpr: 1.0}
  roc_points
end

def round_values_in_curve(curve, *args)
  curve.map{|point|
    point.transform_values{|v| v.round(*args) }
  }.chunk(&:itself).map(&:first) # that's uniq for consequent point
end

def read_motif_hits(stream)
  sequence_lengths = []
  pvalues = []
  stream.each_line.lazy.map(&:chomp).each_slice(2).each{|seq_name, hit_info|
    sequence_length = seq_name.split(':').last.to_i
    pvalue = hit_info.split("\t").first.to_f
    sequence_lengths << sequence_length
    pvalues << pvalue
  }
  {sequence_lengths: sequence_lengths, pvalues: pvalues}
end

show_curve_points = ARGV.delete('--curve-points')
raise 'Specify model length'  unless model_length = Integer(ARGV[0])
raise 'Specify file with model hits or `-` for stdin'  unless control_fn = ARGV[1]

if control_fn == '-'
  motif_hits = read_motif_hits($stdin)
else
  motif_hits = File.open(control_fn){|f| read_motif_hits(f) }
end

pvalues = correct_pvalues(motif_hits[:pvalues],
                          median_sequence_length: median(motif_hits[:sequence_lengths]),
                          model_length: model_length)

roc = roc_curve(pvalues)
logroc = roc
  .reject{|point|
    # we drop all the points which got negative corrected pvalue due to floating point errors
    # and the point (0, 0) because log-transorming of zero is meaningless
    point[:fpr] <= 0
  }
  .map{|point|
    {logfpr: Math.log(point[:fpr]), tpr: point[:tpr]}
  }

roc_auc = calculate_auc(roc, xy_keys: [:fpr, :tpr])
logroc_auc = calculate_auc(logroc, xy_keys: [:logfpr, :tpr])

result = {
  metrics: {
    roc_auc: roc_auc,
    logroc_auc: logroc_auc,
  }
}

result[:supplementary] = {roc: round_values_in_curve(roc, 2), logroc: round_values_in_curve(logroc, 2)}  if show_curve_points

puts result.to_json
