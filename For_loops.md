# For Loops
## Introduction
A for loop is a programming technique that allows us to repeat a command or group of commands across a series of things, such as a list of files, directories, or values.
As a result, they are critical to increasing productivity through automation. 
Using loops, like wildcards and tab completion, minimises the amount of typing necessary (and thus the number of typing errors).
The "for-in" loop is one of the most popular kinds of for loops in Bash.

## Method
Here's a walkthrough of how for loops work in Bash:
```
*for* variable *in* list
do
    # Commands to execute for each item in the list
done
```
The word "for" denotes the beginning of a "For-loop" command.
The term "variable" indicates a variable that will take on the value of each item in the list during each iteration of the loop.
The term "in" separates the variable from the list of things you want to iterate over.
The term "do" denotes the beginning of the job execution list.
The term "done" denotes the completion of a loop.

## Example
```
#!/bin/bash

# Define a list of numbers
numbers=(1 2 3 4 5)

# Iterate over the list and print each number
for num in "${numbers[@]}"
do
    echo "Number: $num"
done
```
In this script:

We build an array called numbers, which has five integer values: 1, 2, 3, 4, and 5.
The for loop iterates across the numbers array.
Inside the loop, the value of the current array member is set to the variable num, and we use echo to output each number with a message.
When you run this script, it will iterate through the list of numbers, printing each one along with the message "Number: " in the terminal.


The shell prompt changes from $ to > and back again as you typing in the loop. 
The second prompt, >, is different to remind you that you haven’t finished typing a complete command yet. 
A semicolon,;, can be used to separate two commands written on a single line.
-
e.g
```
for num in "${numbers[@]}";do echo "Number: $num"; done
```



