import datetime
import sys

import scheduler

date=None

# TODO: Usage: observing_plan.py 2017-11-14
# if no date available: date = today
# Print output filename + other filenames

def validate_date(date_text):
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        raise ValueError("Incorrect data format, should be YYYY-MM-DD.")

if len(sys.argv)==2:
    date=sys.argv[1]
    validate_date(date)

s=scheduler.Scheduler()
s.observing_plan(date=date)
