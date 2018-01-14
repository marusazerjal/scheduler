"""
This is how Jeeves is going to call FunnelWeb Scheduler each time to select the next tile to observe.
"""

import scheduler
reload(scheduler)

# Current telescope position
# This is going to be set by Jeeves
ra_current=None # ideally local meridian
dec_current=-20 # testing

# TODO: Weather input

s=scheduler.Scheduler()
best_tile = s.next_tile(ra_current=ra_current, dec_current=dec_current)
print best_tile

# Should the format of the best_tile be the same as in the observing_plan?

# TODO: if next_tile() fails, read observing_plan.
# Should Scheduler take care of this or Jeeves? Probably Jeeves?
