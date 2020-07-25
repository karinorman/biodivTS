##Pins test
library(tidyverse)
library(pins)

#check github options to make sure files are added to a release and not the repo
getOption("pins.github.release") == 0
#if above expression is FALSE, change .Rprofile option to (usethis::edit_r_profile):
#options(pins.github.release = 0)

#make sure the github PAT is set
#usethis::edit_r_environ()
#GITHUB_PAT="1e03178a12e2c53fb5f640cb30f4473349be3b30"

board_register_github(repo = "karinorman/biodivTS", branch = "master")

#to pin:
#pin(object_name, board = "github")
