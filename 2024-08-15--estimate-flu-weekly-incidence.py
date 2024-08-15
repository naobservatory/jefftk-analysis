#!/usr/bin/env python3
import sys
import csv
import datetime

target_state, target_state_population, start_date, end_date = sys.argv[1:]

target_state_population = int(target_state_population)
parsed_target_start = datetime.date.fromisoformat(start_date)
parsed_target_end = datetime.date.fromisoformat(end_date)

# From p2ra/influenza.py
def parse_mmwr_week(year: int, week: int) -> datetime.date:
    first_day_of_year = datetime.date(year, 1, 1)
    first_sunday_of_year = datetime.date(
        year, 1, 7 - first_day_of_year.weekday()
    )
    assert first_sunday_of_year.weekday() == 6

    if first_day_of_year.weekday() in [3, 4, 5, 6]:
        # First sunday of the year is included in MMWR week 1.
        return first_sunday_of_year + datetime.timedelta(weeks=week - 1)
    else:
        # First sunday of the year is not included in MMWR week 1.
        return first_sunday_of_year + datetime.timedelta(weeks=week - 2)

# https://github.com/naobservatory/p2ra/blob/main/prevalence-data/CDC_WHO_NREVSS_Clinical_Labs.csv
# See https://github.com/naobservatory/p2ra/blob/main/pathogens/influenza.py
# See https://github.com/naobservatory/p2ra/pull/242
def load_cdc_test_data():
    output = {}
    with open("CDC_WHO_NREVSS_Clinical_Labs.csv") as inf:
        cols = None
        for row in csv.reader(inf):
            if row[0].startswith("*"):
                continue  # initial comment

            if not cols:
                cols = row
                continue

            region = row[cols.index("REGION")]
            year = int(row[cols.index("YEAR")])
            mmwr_week = int(row[cols.index("WEEK")])
            parsed_start = parse_mmwr_week(year, mmwr_week)
            centered_date = parsed_start + datetime.timedelta(days=3)
            total_tests = row[cols.index("TOTAL SPECIMENS")]
            positive_a = row[cols.index("TOTAL A")]
            positive_b = row[cols.index("TOTAL B")]

            if total_tests == positive_a == positive_b == "X":
                continue

            if region not in output:
                output[region] = {}

            output[region][centered_date] = (
                int(positive_a),
                int(positive_b),
            )
    return output

cdc_test_data = load_cdc_test_data()

total_positive_tests_a = 0
total_positive_tests_b = 0
total_test_weeks = 0
for centered_date in cdc_test_data[target_state]:
    if parsed_target_start <= centered_date <= parsed_target_end:
        positive_a, positive_b = cdc_test_data[target_state][centered_date]
        total_positive_tests_a += positive_a
        total_positive_tests_b += positive_b
        total_test_weeks += 1

# Now we need to adjust for underreporting

# https://www.cdc.gov/flu-burden/php/data-vis/2023-2024.html
# CDC estimates 35-65M illnesses October 1, 2023, through June 15, 2024
cdc_burden = (35e6 + 65e6) / 2
burden_start_date = datetime.date.fromisoformat("2023-10-01")
burden_end_date = datetime.date.fromisoformat("2024-06-15")
total_a_and_b_within_dates = 0
for state in cdc_test_data:
    for centered_date in cdc_test_data[state]:
        if burden_start_date <= centered_date <= burden_end_date:
            positive_a, positive_b = cdc_test_data[state][centered_date]
            total_a_and_b_within_dates += (positive_a + positive_b)
print("Between %s and %s:" % (burden_start_date, burden_end_date))
print("  CDC counted %s positive tests" % total_a_and_b_within_dates)
print("  CDC estimates %s total infections" % cdc_burden)
underrporting_adjustment = cdc_burden / total_a_and_b_within_dates
print("Each positive test corresponds to an estimated %.1f infections" %
      underrporting_adjustment)
print()
print("Between %s and %s in %s we saw %s positive tests for flu A "
      "and %s for flu B" % (
          start_date, end_date, target_state, total_positive_tests_a,
          total_positive_tests_b))
print("Applying the adjustment for underreporting, there were likely %.0f "
      "flu a infections and %.0f flu b infections." % (
          total_positive_tests_a * underrporting_adjustment,
          total_positive_tests_b * underrporting_adjustment))

print("%s has %s people, so over the course of the season this is %.0f "
      "infections per 100k for flu a, and %.0f for flu b" % (
          target_state, target_state_population,          
     total_positive_tests_a * underrporting_adjustment /
          target_state_population * 100_000,
     total_positive_tests_b * underrporting_adjustment /
          target_state_population * 100_000))
          
print("There were %s weeks between %s and %s, so this is a weekly incidence "
      "per 100k of %.0f for flu a and %.0f for flu b in %s" % (
          total_test_weeks, parsed_target_start, parsed_target_end,
          total_positive_tests_a * underrporting_adjustment /
                    target_state_population * 100_000 / total_test_weeks,
          total_positive_tests_b * underrporting_adjustment /
                    target_state_population * 100_000 / total_test_weeks,
          target_state))
