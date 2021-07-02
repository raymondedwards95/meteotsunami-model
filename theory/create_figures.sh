#!/bin/bash
# Create visualisations for all cases

echo "$(date) - Start creating figures theory" >> last_runs.log

python3 ./theory_part1.py &
python3 ./theory_part2.py &

wait
echo "$(date) - Created figures theory" >> last_runs.log

exit 0
