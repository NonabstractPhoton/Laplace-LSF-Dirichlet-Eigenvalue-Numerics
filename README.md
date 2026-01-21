<!DOCTYPE html>
<html>
  <head>
  </head>
  <body>
      <h1> An Improved Finite Element Based Approach Process for Computing Minimzers to the Dirichlet Eigenvalues </h1>
      <h2> Overview </h2>
      <p>  
      The programs provided in this work are designed to produce computational minimzers for the Dirichlet eigenvalue problem under a volume constraint (|Ω|=1).
      To store information on regions, the evolution algorithm utilizes a dynamic level set approach.
      This provides an advantage over static approaches that require separate rootfinding at each step 
      while still retaining the key benefit of a level set based approach: flexibility in the face of changes in topology.
      The evolution process also implements a unique objective function and corresponding gradient flow that provides resilience 
      in the face of non-simple eigenvalues of arbitrary multiplicity. Finally, an auxillary algorithm was devised to find the optimal
      disjoint union of minimizers to serve as a new minimzer for eigenvalues of greater index.
      </p>
      <h2> User Guide </h2>
      <h4> Candidate</h4>
      <p>
        This is the most important class from the perspective of the end-user, and contains several convenient functions 
        that simplify the inquiry process. The following is a summary of the functions intended for direct use by an end-user.
        More detailed information and source code are available in the file itself.  
        <table>
          <tr>
            <th>Function</th>
            <th>Use</th>
          </tr>
          <tr>
            <td>saveFor</td>
            <td>Saves the given Candidate object as a minimzer for the kth eigenvalue. Optionally provide the output file name.</td>
          </tr>
          <tr>
            <td>view</td>
            <td>Plots the kth minimzer and provides relevant information. Optionally include a specific file to view.</td>
          </tr>
          <tr>
            <td>condense</td>
            <td>Selects a primary minimizer out of Candidates in the folder for the kth eigenvalue.</td>
          </tr>
          <tr>
            <td>clean</td>
            <td>Deletes all but the primary Candidate in the folder for the kth eigenvalue.</td>
          </tr>
          <tr>
            <td>disjointOptimize</td>
            <td>
              Creates and saves a Candidate minimzer for the kth eigenvalue which is the optimal disjoint union 
              of previous minimzers for the minimization of the kth eigenvalue. For best results, primary Candidates should be populated for all
              previous minimizers before running for the kth eigenvalue.
            </td>
          </tr>
          <tr>
            <td>refine</td>
            <td>Starts the evolution process for the kth eigenvalue using its current primary minmizer as initial geometry.</td>
          </tr>
          <tr>
            <td>startFrom</td>
            <td>
              Starts the evolution process for the kth eigenvalue using an input image as initial geometry. 
              This image must grayscaled and use black (0) as the interior of the region.
              This function allows users to draw their own regions and thereby easily test hypotheses.
            </td>
          </tr>
        </table>
      </p>
      <h4>LSFEvolutionDriver</h4>
      <p>
        The driver serves as the launchad from where the user may execute arbitrary jobs. 
        The evolution process and program architecture have been designed to promote parallelization of tasks and thus reduce execution time.
      </p>
      <h4> LaplaceLSF</h4>
      <p>
        This class contains all of the internal code for the dynamic level set used in the evolution process. 
        Instead of the traditional signed distance function, this level set is based on the solution to  ∆u = 1 
        on the region of interest with Dirichlet boundary conditions.
      </p>
      <h4> LaplaceLSFEvolution </h4>
      <p>
        This class is the runtime for the evolution process. Input parameters determine whether the runtime is for 
        manual iteration and closely observing results, or for batch processing of Candidates.
      </p>
  </body>
</html>
