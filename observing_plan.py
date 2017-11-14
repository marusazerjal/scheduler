import datetime
import sys

import scheduler

date=None
time=None

# TODO: Usage: observing_plan.py 2017-11-14 22:45
# if no date and time available: date = today
# Print output filename + other filenames

def validate_date(date_text):
    try:
        datetime.datetime.strptime(date_text, '%Y-%m-%d')
    except ValueError:
        raise ValueError("Incorrect data format, should be YYYY-MM-DD.")

def validate_time(date_text):
    try:
        datetime.datetime.strptime(date_text, '%H:%M')
    except ValueError:
        raise ValueError("Incorrect time format, should be HH:MM.")

if len(sys.argv)==2:
    date=sys.argv[1]
    validate_date(date)
elif len(sys.argv)==3:
    date=sys.argv[1]
    validate_date(date)
    time=sys.argv[2]
    validate_time(time)


s=scheduler.Scheduler()
s.observing_plan(date=date, time=time)
