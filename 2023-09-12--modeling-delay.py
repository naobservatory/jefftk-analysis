#!/usr/bin/env python3

def simulate(doubling_period_weeks, delay_weeks):
  cumulative_infection_threshold = 0.01
  initial_weekly_incidence = 0.000000001
  cumulative_infections = 0
  current_weekly_incidence = 0
  week = 0
  
  delay_weeks += 1

  prev_incidences = [-1]*delay_weeks

  while cumulative_infections < \
        cumulative_infection_threshold:
    week += 1
    current_weekly_incidence = \
        initial_weekly_incidence * 2**(
          week/doubling_period_weeks)
    cumulative_infections += \
        current_weekly_incidence

    prev_incidences.pop(0)
    prev_incidences.append(current_weekly_incidence)

    print(cumulative_infections, current_weekly_incidence)

  #return prev_incidences.pop(0)
  return current_weekly_incidence

#for f in range(50, 500):
for f in [100]:
  doubling_period_weeks = f / 100
  print(doubling_period_weeks,
        simulate(doubling_period_weeks, delay_weeks=0),
        #simulate(doubling_period_weeks, delay_weeks=1),
        #simulate(doubling_period_weeks, delay_weeks=2),
        #simulate(doubling_period_weeks, delay_weeks=4),
        #simulate(doubling_period_weeks, delay_weeks=8),
        sep="\t")

