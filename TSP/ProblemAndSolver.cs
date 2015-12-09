using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;

namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your node data structure and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }


            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        private const int CITY_ICON_SIZE = 5;

        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route  ;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;
        #endregion

        #region Public members
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed); 
            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }


        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        //public void GenerateProblem(int size) // unused
        //{
        //   this.GenerateProblem(size, Modes.Normal);
        //}

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        ///  solve the problem.  This is the entry point for the solver when the run button is clicked
        /// </summary>
        public void branchAndBound()
        {
            //start the clock
            Stopwatch timer = new Stopwatch();
            timer.Start();

            //Find the greatest distance between any two cities
            //O(n^2)
            ArrayList bssfRoute = new ArrayList(); //stores indices of cities in a particular order
            int start = 0;
            int end = 0;
            double maxCost = 0;
            //Find initial vertices to add to subgraph
            for (int k = 0; k < Cities.Length; ++k)
            {
                for (int j = k + 1; j < Cities.Length; ++j)
                {
                    double cost = Cities[k].costToGetTo(Cities[j]);
                    if (maxCost < cost)
                    {
                        start = k;
                        end = j;
                        maxCost = cost;
                    }
                }
            }

            bssfRoute.Add(start);
            bssfRoute.Add(end);
            double minCost;
            int first = 0;
            int last = 0;

            //this will run until are cities have been accounted for
            while (bssfRoute.Count < Cities.Length)
            {
                maxCost = 0;
                start = 0;
                //Find the node that is farthest away from any other node in the sub-graph
                //O(n^2)
                foreach (int cityIndex in bssfRoute)
                {
                    for (int i = 0; i < Cities.Length; ++i)
                    {
                        //Only add if we haven't visited AND if the cost is higher
                        double currentCost = Cities[i].costToGetTo(Cities[cityIndex]);
                        if (!bssfRoute.Contains(i) && currentCost > maxCost)
                        {
                            start = i;
                            maxCost = currentCost;
                        }
                    }
                }

                //Find the edge(i,j) that is the cheapest to insert the newest node into the subgraph
                minCost = double.PositiveInfinity;
                // O(n)
                for (int i = 0; i < bssfRoute.Count; ++i)
                {
                    int i_Index = (int)bssfRoute[i];
                    int j = i + 1;
                    if (j == bssfRoute.Count)
                        j = 0;
                    int j_Index = (int)bssfRoute[j];
                    double cost = findMinimumCost(i_Index, j_Index, start);
                    if (cost < minCost)
                    {
                        first = i_Index;
                        last = j;
                        minCost = cost;
                    }
                }
                //add city to bssf route at the location we just determined
                bssfRoute.Insert(last, start);
            }

            //Convert indexes into Route of cities (nodes)
            int numBssfChanges = 0;
            Route = new ArrayList();
            foreach (int k in bssfRoute) //O(n)
                Route.Add(Cities[k]);

            bssf = new TSPSolution(Route);

            double bssfCost = bssf.costOfRoute();
            bool routeChanged = true;

            //Check to see if swapping any two nodes produces a better (shorter) route 
            while (routeChanged)
            {
                routeChanged = false;
                ArrayList cloneRoute = (ArrayList)Route.Clone();
                //Round robin swap O(nlogn)
                for (int i = 0; i < cloneRoute.Count; ++i)
                {
                    for (int j = i + 1; j < cloneRoute.Count; ++j)
                    {
                        ArrayList tempRoute = (ArrayList)cloneRoute.Clone();
                        City c = (City)tempRoute[i];
                        tempRoute[i] = tempRoute[j];
                        tempRoute[j] = c;
                        TSPSolution newBssf = new TSPSolution(tempRoute);
                        if (bssfCost > newBssf.costOfRoute())
                        {
                            routeChanged = true;
                            cloneRoute = (ArrayList)tempRoute.Clone();
                            bssfCost = newBssf.costOfRoute();
                        }
                    }
                }
                if (routeChanged)
                {
                    Route = (ArrayList)cloneRoute.Clone();
                    bssf = new TSPSolution(Route);
                    numBssfChanges++;
                }
            }
            //update the BSSF to the new route
            bssf = new TSPSolution(Route);
            bssfCost = bssf.costOfRoute();

            //begin preparations for branch and bound
            //create matrix of distances
            double totalCosts = 0;
            double[,] distanceMatrix = new double[Cities.Length, Cities.Length];
            //initialize cost matrix O(n^2)
            for (int i = 0; i < Cities.Length; i++)
            {
                distanceMatrix[i, i] = double.PositiveInfinity;
                for (int j = i + 1; j < Cities.Length; j++)
                {
                    double cost = Cities[i].costToGetTo(Cities[j]);
                    distanceMatrix[i, j] = cost;
                    distanceMatrix[j, i] = cost;
                    totalCosts += cost;
                }
            }

            //get children of start state O(n)
            ArrayList childIndices = new ArrayList();
            for (int i = 1; i < Cities.Length; i++)
                childIndices.Add(i);

            //generate the start state
            State startState = new State(0, childIndices, (double[,])distanceMatrix.Clone(), new ArrayList(), 0, 0, Cities, 0);

            //lower bound on initial state
            double bound = startState.lowerBound;

            //push start state onto priority queue
            PriorityQueue<double, State> myPQ = new PriorityQueue<double, State>();
            myPQ.Enqueue(bound - DEFAULT_SIZE * startState.solutionDepth, startState);

            int maxSize = 0;
            int numStatesCreated = 1;

            //begin branch and bound
            //give B&B a 30 second time limit
            while (!myPQ.IsEmpty && timer.Elapsed.Seconds < 30)
            {
                //keep track of maximun number of states stored at any given iteration
                if (myPQ.Count > maxSize)
                    maxSize = myPQ.Count;

                State u = myPQ.Dequeue().Value;

                //lazy pruning
                if (u.lowerBound < bssfCost)
                {
                    ArrayList childStates = u.generateChildStates();
                    numStatesCreated += childStates.Count;
                    foreach (State s in childStates)
                    {
                        if (timer.Elapsed.Seconds > 30)
                            break;
                        //only deal with child state if bound is lower than bssf
                        if (s.lowerBound < bssfCost)
                        {
                            if (s.solutionFound && s.currentCost < bssfCost)
                            {
                                bssfCost = s.currentCost;
                                bssfRoute = s.pathSoFar;
                            }
                            else
                            {
                                //double cBound = s.lowerBound;
                                myPQ.Enqueue(s.lowerBound - DEFAULT_SIZE * s.solutionDepth, s);
                            }
                            //numBssfChanges++;
                        }
                    }
                }
            }
            timer.Stop();
            bssf = new TSPSolution(Route);
            //bssf = apply2Opt(bssf);

            //Console.WriteLine("distance from 1 to 6: {0}", ((City)bssf.Route[1]).costToGetTo((City)bssf.Route[6]));
            //Console.WriteLine("distance from 2 to 7: {0}", ((City)bssf.Route[2]).costToGetTo((City)bssf.Route[7]));

            Console.WriteLine("Max number of states stored = {0}", maxSize);
            Console.WriteLine("bssf changed {0} times", numBssfChanges);
            Console.WriteLine("States created = {0}", numStatesCreated);
            Console.WriteLine("Number of states pruned = {0}", (numStatesCreated - maxSize));
            
            //write out solution
            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            TimeSpan totalTime = timer.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}", totalTime.Minutes, totalTime.Seconds);
            Program.MainForm.tbElapsedTime.Text = elapsedTime;
            Program.MainForm.Invalidate();

            Console.WriteLine(ToStringSolution());
        }
        public void TwoOpt()
        {
            //start the clock
            Stopwatch timer = new Stopwatch();
            timer.Start();

            //Find the greatest distance between any two cities
            //O(n^2)
            ArrayList bssfRoute = new ArrayList(); //stores indices of cities in a particular order
            int start = 0;
            int end = 0;
            double maxCost = 0;
            //Find initial vertices to add to subgraph
            for (int k = 0; k < Cities.Length; ++k)
            {
                for (int j = k + 1; j < Cities.Length; ++j)
                {
                    double cost = Cities[k].costToGetTo(Cities[j]);
                    if (maxCost < cost)
                    {
                        start = k;
                        end = j;
                        maxCost = cost;
                    }
                }
            }

            bssfRoute.Add(start);
            bssfRoute.Add(end);
            double minCost;
            int first = 0;
            int last = 0;

            //this will run until are cities have been accounted for
            while (bssfRoute.Count < Cities.Length)
            {
                maxCost = 0;
                start = 0;
                //Find the node that is farthest away from any other node in the sub-graph
                //O(n^2)
                foreach (int cityIndex in bssfRoute)
                {
                    for (int i = 0; i < Cities.Length; ++i)
                    {
                        //Only add if we haven't visited AND if the cost is higher
                        double currentCost = Cities[i].costToGetTo(Cities[cityIndex]);
                        if (!bssfRoute.Contains(i) && currentCost > maxCost)
                        {
                            start = i;
                            maxCost = currentCost;
                        }
                    }
                }

                //Find the edge(i,j) that is the cheapest to insert the newest node into the subgraph
                minCost = double.PositiveInfinity;
                // O(n)
                for (int i = 0; i < bssfRoute.Count; ++i)
                {
                    int i_Index = (int)bssfRoute[i];
                    int j = i + 1;
                    if (j == bssfRoute.Count)
                        j = 0;
                    int j_Index = (int)bssfRoute[j];
                    double cost = findMinimumCost(i_Index, j_Index, start);
                    if (cost < minCost)
                    {
                        first = i_Index;
                        last = j;
                        minCost = cost;
                    }
                }
                //add city to bssf route at the location we just determined
                bssfRoute.Insert(last, start);
            }

            //Convert indexes into Route of cities (nodes)
            int numBssfChanges = 0;
            Route = new ArrayList();
            foreach (int k in bssfRoute) //O(n)
                Route.Add(Cities[k]);

            bssf = new TSPSolution(Route);

            double bssfCost = bssf.costOfRoute();
            bool routeChanged = true;

            //Check to see if swapping any two nodes produces a better (shorter) route 
            while (routeChanged)
            {
                routeChanged = false;
                ArrayList cloneRoute = (ArrayList)Route.Clone();
                //Round robin swap O(nlogn)
                for (int i = 0; i < cloneRoute.Count; ++i)
                {
                    for (int j = i + 1; j < cloneRoute.Count; ++j)
                    {
                        ArrayList tempRoute = (ArrayList)cloneRoute.Clone();
                        City c = (City)tempRoute[i];
                        tempRoute[i] = tempRoute[j];
                        tempRoute[j] = c;
                        TSPSolution newBssf = new TSPSolution(tempRoute);
                        if (bssfCost > newBssf.costOfRoute())
                        {
                            routeChanged = true;
                            cloneRoute = (ArrayList)tempRoute.Clone();
                            bssfCost = newBssf.costOfRoute();
                        }
                    }
                }
                if (routeChanged)
                {
                    Route = (ArrayList)cloneRoute.Clone();
                    bssf = new TSPSolution(Route);
                    numBssfChanges++;
                }
            }
            //update the BSSF to the new route
            bssf = new TSPSolution(Route);
            bssfCost = bssf.costOfRoute();

            //begin preparations for branch and bound
            //create matrix of distances
            double totalCosts = 0;
            double[,] distanceMatrix = new double[Cities.Length, Cities.Length];
            //initialize cost matrix O(n^2)
            for (int i = 0; i < Cities.Length; i++)
            {
                distanceMatrix[i, i] = double.PositiveInfinity;
                for (int j = i + 1; j < Cities.Length; j++)
                {
                    double cost = Cities[i].costToGetTo(Cities[j]);
                    distanceMatrix[i, j] = cost;
                    distanceMatrix[j, i] = cost;
                    totalCosts += cost;
                }
            }

            //get children of start state O(n)
            ArrayList childIndices = new ArrayList();
            for (int i = 1; i < Cities.Length; i++)
                childIndices.Add(i);

            //generate the start state
            State startState = new State(0, childIndices, (double[,])distanceMatrix.Clone(), new ArrayList(), 0, 0, Cities, 0);

            //lower bound on initial state
            double bound = startState.lowerBound;

            //push start state onto priority queue
            PriorityQueue<double, State> myPQ = new PriorityQueue<double, State>();
            myPQ.Enqueue(bound - DEFAULT_SIZE * startState.solutionDepth, startState);

            int maxSize = 0;
            int numStatesCreated = 1;

            //begin branch and bound
            //give B&B a 30 second time limit
            while (!myPQ.IsEmpty && timer.Elapsed.Seconds < 30)
            {
                //keep track of maximun number of states stored at any given iteration
                if (myPQ.Count > maxSize)
                    maxSize = myPQ.Count;

                State u = myPQ.Dequeue().Value;

                //lazy pruning
                if (u.lowerBound < bssfCost)
                {
                    ArrayList childStates = u.generateChildStates();
                    numStatesCreated += childStates.Count;
                    foreach (State s in childStates)
                    {
                        if (timer.Elapsed.Seconds > 30)
                            break;
                        //only deal with child state if bound is lower than bssf
                        if (s.lowerBound < bssfCost)
                        {
                            if (s.solutionFound && s.currentCost < bssfCost)
                            {
                                bssfCost = s.currentCost;
                                bssfRoute = s.pathSoFar;
                            }
                            else
                            {
                                //double cBound = s.lowerBound;
                                myPQ.Enqueue(s.lowerBound - DEFAULT_SIZE * s.solutionDepth, s);
                            }
                            //numBssfChanges++;
                        }
                    }
                }
            }
            timer.Stop();
            bssf = new TSPSolution(Route);
            bssf = apply2Opt(bssf);

            //Console.WriteLine("distance from 1 to 6: {0}", ((City)bssf.Route[1]).costToGetTo((City)bssf.Route[6]));
            //Console.WriteLine("distance from 2 to 7: {0}", ((City)bssf.Route[2]).costToGetTo((City)bssf.Route[7]));

            Console.WriteLine("Max number of states stored = {0}", maxSize);
            Console.WriteLine("bssf changed {0} times", numBssfChanges);
            Console.WriteLine("States created = {0}", numStatesCreated);
            Console.WriteLine("Number of states pruned = {0}", (numStatesCreated - maxSize));

            //write out solution
            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            TimeSpan totalTime = timer.Elapsed;
            string elapsedTime = string.Format("{0:00}:{1:00}", totalTime.Minutes, totalTime.Seconds);
            Program.MainForm.tbElapsedTime.Text = elapsedTime;
            Program.MainForm.Invalidate();

            Console.WriteLine(ToStringSolution());
        }

        public double findMinimumCost(int i, int j, int k)
        {
            double cost = Cities[k].costToGetTo(Cities[i]) +
                Cities[j].costToGetTo(Cities[k]);
            cost -= Cities[i].costToGetTo(Cities[j]);
            return cost;
        }

        public string ToStringSolution()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append('[');

            foreach (City city in Route)
            {
                sb.Append("[" + city.X + "," + city.Y + "],");
            }
            //Remove last comma
            sb.Remove(sb.Length - 1, 1); 
            sb.Append(']');

            return sb.ToString();
        }

        private TSPSolution apply2Opt(TSPSolution myRoute)
        {
            if (myRoute.Route.Count < 4)
                return myRoute;

            //double[,] distanceMatrix = new double[myRoute.Route.Count, myRoute.Route.Count];
            ////initialize cost matrix O(n^2)
            //for (int i = 0; i < myRoute.Route.Count; i++)
            //{
            //    distanceMatrix[i, i] = double.PositiveInfinity;
            //    for (int j = i + 1; j < myRoute.Route.Count; j++)
            //    {
            //        double cost = ((City)myRoute.Route[i]).costToGetTo((City)myRoute.Route[j]);
            //        distanceMatrix[i, j] = cost;
            //        distanceMatrix[j, i] = cost;
            //        //totalCosts += cost;
            //    }
            //}

            double minChange = -1;
            while (minChange < 0)
            {
                minChange = 0;
                int min_i = 0;
                int min_j = 0;
                for (int i = 0; i < myRoute.Route.Count - 2; i++)
                {
                    for (int j = i+2; j < myRoute.Route.Count - 1; j++)
                    {
                        double change = ((City)myRoute.Route[i]).costToGetTo((City)myRoute.Route[j]);
                        change += ((City)myRoute.Route[i + 1]).costToGetTo((City)myRoute.Route[j + 1]);
                        change -= ((City)myRoute.Route[i]).costToGetTo((City)myRoute.Route[i + 1]);
                        change -= ((City)myRoute.Route[j]).costToGetTo((City)myRoute.Route[j + 1]);


                        if (minChange > change && isRouteShorter(change, myRoute, i, j))
                        {
                            minChange = change;
                            min_i = i;
                            min_j = j;
                        }
                    }
                }
                // apply min_i/min_j move
                bool middleCanBeReversed = true;
                for (int i = min_j; i > min_i; i--)
                    if (((City)myRoute.Route[i]).costToGetTo((City)myRoute.Route[i - 1]).Equals(double.PositiveInfinity))
                        middleCanBeReversed = false;
                if (middleCanBeReversed)
                    myRoute = twoOptSwap(myRoute, min_i, min_j);
            }
            return myRoute;
        }

        private bool isRouteShorter(double change, TSPSolution route, int i, int j)
        {
            double originalLength = 0;
            double newLength = 0;
            for (int k = i; k < j; k++)
                originalLength += ((City)route.Route[k]).costToGetTo((City)route.Route[k + 1]);
            for (int k = j; k > i; k--)
                newLength += ((City)route.Route[k]).costToGetTo((City)route.Route[k - 1]);
            newLength += change;

            return (newLength < originalLength);
        }

        private TSPSolution twoOptSwap(TSPSolution route, int i, int j)
        {
            //function 2optSwap(route, i, j)
            //Input: TSP solution route; nodes i and j in the route, 1 < i < j
            //Output: A new TSP route with the swap

            TSPSolution new_route = new TSPSolution(new ArrayList());
            for (int k = 0; k < i + 1; k++)// k = 1 to i - 1)
            {
                new_route.Route.Add(route.Route[k]);
                //add(new_route, route[k])
            }
            for (int k = j; k > i; k--)// k = j to i:
            {
                new_route.Route.Add(route.Route[k]);
                //add(new_route, route[k])
            }
            for (int k = j + 1; k < route.Route.Count; k++) //k = j + 1 to size(route):
            {
                new_route.Route.Add(route.Route[k]);
                //add(new_route, route[k])
            }

            return new_route;

        }

        public void solveProblemGreedy()
        {

            Route = new ArrayList();
            Route.Add(Cities[0]);
            while (Route.Count < Cities.Length)
            {  //TODO do hard mode
                Route.Add(getClosestTo((City)Route[Route.Count - 1], Route));
            }
            // call this the best solution so far.  bssf is the route that will be drawn by the Draw method. 
            bssf = new TSPSolution(Route);
            // update the cost of the tour. 
            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            // do a refresh. 
            Program.MainForm.Invalidate();
        }

        public City getClosestTo(City c, ArrayList routeSoFar)
        {
            City minCity = null;
            double minDist = Double.MaxValue;
            for (int x = 0; x < Cities.Length; x++)
            {
                if (c.costToGetTo(Cities[x]) < minDist && Cities[x] != c && !routeSoFar.Contains(Cities[x]))
                {
                    minCity = Cities[x];
                    minDist = c.costToGetTo(Cities[x]);
                }
            }
            if (minDist == Double.PositiveInfinity)
            {
                throw new badPathException();
            }

            return minCity;

        }

        private class badPathException : Exception
        {
            public badPathException()
            {
            }

            public badPathException(string message)
                : base(message)
            {
            }

            public badPathException(string message, Exception inner)
                : base(message, inner)
            {
            }
        }

        public void solveRandomPath()
        {
            Route = new ArrayList();
            //Route.Capacity = Cities.Length;
            Random rand = new Random();
            // this is the trivial solution. 
            for (int x = 0; x < Cities.Length; x++)
            {
                Route.Add(Cities[x]);
            }

            //fisher-yates shuffle
            for (int i = Route.Count - 1; i > 0; i--)
            {
                int index = rand.Next(i + 1);
                // Simple swap
                City c = (City)Route[index];
                Route[index] = Route[i];
                Route[i] = c;
            }

            // call this the best solution so far.  bssf is the route that will be drawn by the Draw method. 
            bssf = new TSPSolution(Route);
            // update the cost of the tour. 
            Program.MainForm.tbCostOfTour.Text = " " + bssf.costOfRoute();
            // do a refresh. 
            Program.MainForm.Invalidate();

        }

        #endregion
    }

}
