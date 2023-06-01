# Least-Cost path analysis of camel caravan trade routes in ancient Arabia
A geospatial analysis used in my Masterthesis. Its a Least Cost Path analysis based around the R package gdistance. It employs different cost algorithms (Tobler 1994, Minetti 2002, Herzog 2012) to calculate cost surfaces and provide comparability between them. On them Dijkstra´s (1959) pathfinding algorithm is applied to identify the Least Cost path. Also utilises the gdistance function to generate passage raster to explore movement corridors, a better approximation to human movement than linear features

## Literature

the research is published as an article (https://brill.com/display/book/9789004527119/BP000015.xml). A publication as a monography is in progress

## Tech stack

written in R, used Rstudio as IDE, used the gdistance library for the actual LCP functions

## 5 years later

funny to look at it 5 years later. code-wise it is pretty rudimentary. back then i was not able to write loops so i just chained the code for all of 20ish routes. biggest problem back then was that all the operation run in the RAM and for 200MB DEM's that accumulated to 70GB for the transition matrix. rare to find a PC specced like this but i was lucky. Today i would probably look for a library that has a different memory management at the cost of increased calculation time. or consider hosted HPC services
