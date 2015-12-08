using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TSP
{
    class State
    {
        //Definition of state
        public ArrayList pathSoFar;
        public int currIndex;
        public ArrayList childrenIndices;
        public double[,] costMatrix;
        public double currentCost;
        public City[] cities;
        public double lowerBound;
        public int solutionDepth;
        public bool solutionFound;

        public State(double lowerBound, ArrayList childIndices, double[,] costMatrix, ArrayList pathSoFar,
            double startCost, int currIndex, City[] cities, int depth)
        {
            this.lowerBound = lowerBound;
            this.childrenIndices = childIndices;
            this.costMatrix = costMatrix;
            this.pathSoFar = pathSoFar;
            this.solutionFound = false;
            this.currentCost = startCost;
            this.currIndex = currIndex;
            this.cities = cities;
            this.solutionDepth = depth;

            //a solution is found when there are no more children to explore 
            if (childrenIndices.Count == 0)
                solutionFound = true;

            //generate the lower bound
            this.lowerBound = lowerBound + reducedCostMatrix();

            if (solutionFound)
            {
                currentCost += cities[0].costToGetTo(cities[currIndex]);
                this.pathSoFar.Add(currIndex);
            }
        }

        public ArrayList generateChildStates()
        {
            ArrayList children = new ArrayList();
            //create child state for each possibility
            foreach (int i in childrenIndices)
            {

                ArrayList newChildren = (ArrayList)childrenIndices.Clone();
                newChildren.Remove(i);

                // update the path based on the current child
                ArrayList newPath = (ArrayList)pathSoFar.Clone();
                newPath.Add(currIndex);
                double cost = cities[currIndex].costToGetTo(cities[i]);

                //reduced matrix for new state
                double[,] newCost = (double[,])costMatrix.Clone();
                for (int j = 0; j <= newCost.GetUpperBound(0); j++)
                    newCost[j, currIndex] = double.PositiveInfinity;

                children.Add(new State(lowerBound + costMatrix[currIndex, i],
                    newChildren, newCost, newPath,
                    currentCost + cost, i, cities, solutionDepth + 1));
            }
            return children;
        }

        private double reducedCostMatrix()
        {
            double[,] reducedMatrix = costMatrix;
            double costSoFar = 0;
            double rowsNotChanged = 0;

            //row reduction
            ArrayList rows = (ArrayList)childrenIndices.Clone();
            rows.Add(currIndex);
            foreach (int i in rows)
            {
                double minCost = double.PositiveInfinity;
                foreach (int j in childrenIndices)
                    if (reducedMatrix[i, j] < minCost)
                        minCost = reducedMatrix[i, j];

                if (!double.IsPositiveInfinity(minCost))
                {
                    foreach (int j in childrenIndices)
                        reducedMatrix[i, j] -= minCost;
                    costSoFar += minCost;
                }
                else
                    rowsNotChanged++;
            }

            //column reduction
            foreach (int i in childrenIndices)
            {
                double minCost = double.PositiveInfinity;
                foreach (int j in rows)
                    if (reducedMatrix[j, i] < minCost)
                        minCost = reducedMatrix[j, i];
                if (!double.IsPositiveInfinity(minCost))
                    foreach (int j in rows)
                        reducedMatrix[j, i] -= minCost;
                else
                    rowsNotChanged++;
            }

            if (rowsNotChanged >= reducedMatrix.GetUpperBound(0))
            {
                solutionFound = true;
                costSoFar += cities[1].costToGetTo(cities[currIndex]);
            }

            if (solutionFound)
                return 0;

            return costSoFar;
        }

    }
}
