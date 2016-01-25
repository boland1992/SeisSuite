The following is a list of things to-do to improve the SeisSuite package.

1. Attempt to keep the CHANGES.txt file up-to-date
2. Allow for choice of where the config is stored in case the os.getcwd() option is not desirable. 
3. Allow for choice of different extensions in create_database.py
4. Get the new_station_search.py script to work with config file!
5. Integrate option to either have the linear stack or the phase-weighted stack used for FTAN
6. Add if __name__ == '__main__' support for all tools that don't already have it. 
7. Improve the initialisation of the SQL databases from 01_database_init.py to run as functions and not run during importation. 
8. Create functional tool to search frequency response window in the response.db, and also create tests to check if the response.db is correct!
9. Create new timeline database to keep track of when the processing for preprocessing and xcorr is. This is in case the whole terminal crashes or is interupted for any reason. Currently it uses pickle, but SQL feels like a better fit for the task. 
10. Fix problems with False multiprocess option in 02_timeseries_process.py
11. Fix problem with not being able to set either nothing or / at the end of the FOLDER path in .cnf
12. Create SNR with time for all station pairs' signals and for each stack type e.g. SNR_lin stack or SNR_pws (phase-weighted Stack). 
13. Create tool for showing SNR calculation technique. Visualise the 'noise' vs. 'signal' sections of the waveform.
14. Allow for multiple runs with varyious different config files!
15. 
16. Impliment total processing time.
17. Automate a text output that states total processing time, number of stations, number of station pairs, input configuration file paremeters (including which preprocessing and post-processing methods were used) and the related SNR information.  
18. Create functions for all normalisation types and options to use each of them
19. Create functions for all SNR types, find a way to test which is best? and reference them all!
20. Set alternative for downsampling to max as min sample rate for database! currently the application
will NOT work if DOWNSAMPLE = False in config file.
21. Allow preprocessing steps to be set in an individual order (i.e. function construction!)
22. rename config check function to something a little more appropriate and allow for user input of both lists and individual settings.
23. Try creating a NEW config text/ anything file to save current parsed configuration global information in!
24. have a way of accessing the current configuration file location! or something like that. 
25. allow for automatic SQL database file structuring OR alternatives e.g. BUDS
26. Change output print for when xcorr sets are saving. Give option for individual savings 
27. Change re-start pickle to output to xcorr
28. changing the line here to show up on new branch


Fixed or new:

1. Fixed issues with pspreprocess not yielding correct results!
2. Tidied up the SNR_table function in pscrosscorr, it is now working nominally. 
3. Managed to save incremental time-steps in (2,N) shaped array for the SNR vs. time graphs. 
4. Preprocessing steps can now be set True or False in config files!
