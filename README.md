# Monte-Carlo
Code for weiner process and Monte Carlo simulations with Poisson jumps, csv dump, csv readfile, ARCH volatility adjustments, and plotting

To use, you will have to manually input the appropriate stock information. Make sure that you either turn off the poisson jumps or enter the appropriate values. When you input your volatility spreadsheet, try to format it the same as the one I have uploaded, or just import it yourself. 

Model should run reasonably fast, remember that a good model should have at least 500 time steps and 2000 trials. This can take a few minutes to run. Anything over 5 minutes and you may want to make sure you aren't getting stuck in a loop.

The finite difference call pricing model is there to provide another way to value a call option. Using similar inputs, you should get similar results. The monte carlo model will be more customizable to add in the randomness of the markets, jumps in prices, and improved accuracy of volatility. Again, you just have to enter the values for the variables and run it in a python console.
