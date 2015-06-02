

# Calculate the similarity using the default method of optimal string alignment
stringsim("ca", "abc")

# Calculate the similarity using the Jaro-Winkler method
# The p argument is passed on to stringdist
stringsim('MARTHA','MATHRA',method='jw', p=0.1)

