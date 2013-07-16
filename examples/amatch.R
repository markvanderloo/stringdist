
# lets see which sci-fi heroes are stringdistantly nearest
amatch("leia",c("uhura","leela"))

# we can restrict the search
amatch("leia",c("uhura","leela"),maxDist=1)

# this is always true (maxDist is Inf by default)
ain("leia",c("uhura","leela"))

# so we set a reasonable maxDist
ain("leia",c("uhura","leela"), maxDist=1)


