# LeastSquaresIntroToML
Project introducing ML math concepts via least-squares curve fitting

Hello,

Last year, I noticed that some people at a local AI Meetup appeared to be new to ML and seemed to have trouble following some of the more advanced presentations (as did I!).  Some of these people may have had rusty Calculus knowledge as well (me too!).

So, I decided to create a small project introducing them to least-squares curve fitting (primarily cost functions, models, Gradient Descent) and also to refresh their memory of some relevant concepts from Calculus (derivative, slope, gradient).  Since curve fitting has concepts similar to those in ML and because the math formulation can be presented with less formality, I felt this was a good/gentle introduction to ML concepts.

I decided to use R for this project since I had taken a beginner class in R in the Johns Hopkins Coursera Data Science sequence several years ago, and seemed to recall that R's plotting and numerical computations similar to cost functions were simple, mostly because R is in inherently 'vectorized.'  So in doing this project, I had to refresh my minimal prior R knowledge from years ago, and also learn its model-fitting and plotting APIs.

Note: I am *not* an R expert, so the included code is as-is.  The code was intended to make the computations concrete.

In order to make the code more readable to a wider audience of potentially non-R coders, I chose to limit the use of R idioms and also used the '=' assignment operator instead of the '<-' operator.  In addition, I used camel-casing instead of dashes in variable and function names.

In addition to this repo on Github, I have a web page for this project on my personal website(URL at bottom) which has more information about the code/R code techniques used, and a grid/table of many of the plot thumbnails which allows the user to click on a plot thumbnail and be taken to a sub-page with both the R code and the associated full-sized plot PNG file.

The code consists of (2) main R files:
  - Examples.R [example-specific computations and plotting code]
  - Util.R     [a basic library of helper functions I developed]

In addition:
  - Each individual plot has its own sub-directory containing a plot PNG file and an R file.
  - The per-plot R files contain only the code necessary for that plot (self-contained).

Finally, I am including the PowerPoint slide deck .pptx file and a PDF export of this file.  In addition, I split the deck into (2) parts, and these split deck files are also included.

Note: This work is copyrighted.  Contact me via email (below) for usage inquiries.

Thanks for looking.  Feedback is welcome via my Gmail/email account, below.

Richard Creamer, 1/20/2020
Email: 2to32minus1@gmail.com
Website: https://sites.google.com/site/rickcreamer
Project page: https://sites.google.com/site/rickcreamer/richard-creamers-personal-website/ml-curve-fitting-intro
