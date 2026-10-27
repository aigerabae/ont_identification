Test run:
```
cd ~/biostar/NCB/hp_august2
conda env create -f environment.yml && conda activate plant_its
bash tests/run_test.sh        # optional: about 5 minutes, no internet; should print the test table
```

Building database (once):
```
cd ~/biostar/NCB/hp_august2
conda activate plant_its
chmod +x db/*.sh db/*.py scripts/*.sh

bash tests/run_test.sh          # optional sanity check, ~5 min; wait for the final table
bash db/get_taxdump.sh          # must finish first: creates db/taxdump/nodes.dmp
bash db/build_custom_db.sh      # hours; run it in tmux or screen so it survives a closed terminal
```
