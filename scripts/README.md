### Generate words 
1. Create a directory where the output files will go.
2. Copy the `global.sec.word` file for each lecture (right now it's only for 18.06-1999/L02) from the exps folder 
3. To generate the mapping between plus and words, run ./word_discovery.sh

### Warning
In `word_discovery.sh`, a binary called `ctl_exec` is called. This is a tool the SLS group at CSAIL MIT uses to run jobs in parallel. You have to replace this or just modify the `plu_to_wrd_in_samples.ctl` file so that you can run them locally.
