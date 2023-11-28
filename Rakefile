require 'bundler/setup'
$:.unshift File.expand_path('../lib', __FILE__)

# rake spec
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) { |t| t.verbose = false   }

# rake console
task :console do
  require 'pry'
  ARGV.clear
  Pry.start
end

## data/simdata contains  Hc, Ha, 18xHn datasets for one region
## data/simdata_v2 contains  Hc, Ha, Hc2, Ha2, 18xHn datasets for one region

desc "simulate Hc, Ha, Hc2, Ha2, 18xHn datasets for one region. x 20 in each run. x 1000 runs = 20000 datasets"
task :sim do
  system("qR.rb -r -y 0-999 sim_data.R --args d=data/simdata_v2")
end

## final combined data file
collate_file="data/simdata_v2/collated.RData"

desc "collate"
task :collate => collate_file
file collate_file do
  system("qR.rb -r collate.R --args d=data/simdata_v2 outfile=#{collate_file}")
end
