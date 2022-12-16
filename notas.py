#!/usr/bin/env python3

import re

# Set variables to keep count on the sum of grades ("total")  and the number of inserted grades ("grade_count")
total = 0
grade_count = 0

# Start a while loop to get grades from the user
while grade_count < 10:
    grade = input("\nEnter a grade:")

    # Check whether the user input consists of floating number. If not, print an error message and go back to the start of the while loop
    try:
        grade = float(grade)
    except:
        print("Only numeric characters are accepted. Please, try again.")
        continue

    # If the user entry is valid, sum it to the total, update the grade count and print grade_count so the user can keep track on the number of inputted grades
    total += float(grade)
    grade_count += 1
    print("Grade count:", grade_count)

# Assign the average to average_grade
average_grade = round(float(total/10), 2)

# Print the average grade
print("\nAverage grade:", average_grade)
