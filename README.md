# Compile 
A basic makefile is provided as a starting point for compiling the code. You will want to edit the makefile (for example, specify the correct path of the compiler) to compile the code.

# Source Code
Please note that a large portion of the source code in this repo was contributed by Mark Johnson (http://web.science.mq.edu.au/~mjohnson/Software.htm). I'm very grateful for that he made his source code for his NAACL'09 paper avaialble, and please cite his work if you find his portion of the code helpful for your own project.

# File Description
For files that were contributed by me, I'm providing a brief description for each of them here. If you find the description unclear or have any questions, please feel free to send me an email at chiaying@csail.mit.edu or file a pull request to change it!

To better understand the descriptions, I suggest the reader to refer to this paper (http://people.csail.mit.edu/chiaying/publications/acl2012.pdf). I'll use terminologies that are consistent with those used in that paper.

* bound.cc: represents the boundary variables. A bound object consists of a few speech feature frames. 
* segment.cc: represents the segments. A segment object consists of one or more bound objects.
* cluster.cc: represents an automatically discovered unit. You can think of each cluster as a phone unit. 
* config.cc: sets up all the experiment configuration.
* datum.cc: represents an input speech utterance.
* gmm.cc: an implementation of a Gaussian mixture model.
* mixture.cc: a Gaussian mixture component.
* sampler.cc: contains code for doing the sampling-based inference steps.
