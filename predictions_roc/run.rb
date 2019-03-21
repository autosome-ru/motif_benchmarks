require 'json'

def calculate_auc(points, xy_keys: [0, 1])
  xy_coords = points.map{|point| point.values_at(*xy_keys) }
  xy_coords.sort_by{|x,y|
    x
  }.each_cons(2).map{|(x1, y1), (x2, y2)|
    (x2 - x1) * (y1 + y2) / 2.0
  }.inject(0.0, &:+)
end

def roc_curve(ordered_labels)
  num_positive = ordered_labels.count('+')
  num_negative = ordered_labels.count('-')

  tp, fp = 0, 0
  roc_points = []
  roc_points << {tpr: 0.0, fpr: 0.0}
  ordered_labels.each_with_index{|label, index|
    if label == '+'
      tp += 1
    else
      fp += 1
    end
    roc_points << {tpr: tp.to_f / num_positive, fpr: fp.to_f / num_negative}
  }
  roc_points << {tpr: 1.0, fpr: 1.0}
  roc_points
end

def round_values_in_curve(curve, *args)
  curve.map{|point|
    point.transform_values{|v| v.round(*args) }
  }.chunk(&:itself).map(&:first) # that's uniq for consequent point
end

config_fn = "/workdir/config.json"
result_fn = "/workdir/persistent/result.json"
config = JSON.parse(File.read(config_fn))

ground_truth_fn = '/benchmark_specific_data/ground_truth.txt'
labels = File.readlines(ground_truth_fn).map(&:strip)

predictions = config['predictions']
ordered_labels = labels.zip(predictions).sort_by{|label, prediction| -prediction }.map{|label, prediction| label }

roc = roc_curve(ordered_labels)
roc_auc = calculate_auc(roc, xy_keys: [:fpr, :tpr])

results = {
  metrics: {
    roc_auc: roc_auc,
  },
  supplementary: {
    roc: round_values_in_curve(roc, 2),
  },
}

File.write(result_fn, results.to_json)
