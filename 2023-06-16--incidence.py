#!/usr/bin/env python3

#daily_growth_rates = [0.05, 0.07, 0.1, 0.15, 0.2, 0.25]

SUCCEPTIBLE=-1
RECOVERED=-2

population = [SUCCEPTIBLE] * 300
population[0] = 0 # one person infected on day zero


# Each person is infectious for this many days.
INFECTIOUS_PERIOD_DAYS=7
# Each person infects this many other people daily, unless they've already been
# infected.
DAILY_INFECTIONS = 0.5

day = 0
while True:
    day += 1
    for i in range(population):
        if population[i] >= 0:
            if day - population[i] > INFECTIOUS_PERIOD_DAYS:
                population[i] = RECOVERED
            else:
                

while infected + recovered < 0.01 * population:
    infected *= 2


        
        pct_infected *= (1 + daily_growth_rate)
        days += 1
    print("%.0f%% to have infected 1%% of people: %s days" % (
        daily_growth_rate * 100, days))
