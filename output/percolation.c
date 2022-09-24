
/* ------------------------------------------------------------------------ */
/*        Site-Percolation Model on the Simple 2D Lattice                   */
/*        Free Boundary Conditions                                          */
/*        Cluster Analysis with the Hoshen-Koppelman-Algorithm              */         
/*                                                                          */
/*        (c) Dieter W. Heermann                                            */
/*                                                                          */
/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
/*                                                                          */
/*        Required parameters for the simulation are:                       */
/*        -------------------------------------------                       */
/*                                                                          */
/*        l_size    linear system size                                      */ 
/*        pB        percolation probability                                 */
/*        mySeed    seed for the random number generator.                   */
/*        mcsmax    maximum number of monte carlo steps                     */
/*                  (sampling size)                                         */
/* ------------------------------------------------------------------------ */

/* ==== Program parameters ==== */

# define l_size          128
# define pB              0.6
# define sampling_size   100
# define mySeed            4711

/* ==== Include files ==== */

# include       <math.h>
# include       <stdio.h>

int init_r250(int seed, int *m_f_ptr);
int r250 (int n, float *x_ptr, int *m_f_ptr,int save);

int main()
  {
    /* -------------------------------------------------------------------- */
    /*                        Declarations                                  */
    /* -------------------------------------------------------------------- */

    /*   ==== Arrays    ==== */

    int iss[l_size][l_size];       /* 2D site percolation lattice          */
    int row[l_size+1];             /* temp vector                          */
    int lptr[l_size * l_size];     /* Contains the pointers. The dimension */
                                   /* is set, such that a concentration of */
                                   /* one, i.e., one huge cluster is       */
                                   /* possible.                            */
    int mf[251];                   /* contains numbers of the RNG          */

    float vec1[l_size * l_size];
    float nsc[l_size * l_size];
    float ran[l_size+1];           /* holds the random numbers             */

    int nofs, irow, i;
    int la, cl, ms, ix, iy, up, mlarge, mcsmax;
    int iseed, large;
    int mcs;
    int icol, left, mini, maxi;
    int nis;
    int l_sqr;
    int hold;
    int return_code;
    int    c1,c2;
    float r,d;
    int   *mf_ptr;
    float *ran_ptr;
    float pAverage;
    float ri, rj;
    int largeLabel;

    /* -------------------------------------------------------------------- */
    /*                       End of Declarations                            */
    /* -------------------------------------------------------------------- */

    mcsmax = sampling_size;
    iseed = mySeed;

    /* -------------------------------------------------------------------- */
    /*                    Set and Initialize                                */
    /* -------------------------------------------------------------------- */

 
    r = 0.01;                  /* radius of circles, discs, ..              */ 
    d = 1.5 * r;               /* width and height of rectangles            */

    /*  ==== set up the random number generator ==== */

    iseed       = (iseed << 2) + 1;
    mf_ptr      = mf;
    ran_ptr     = ran;
    hold        = 1;
    return_code = init_r250(iseed,mf_ptr);

    /* ==== Initialize the label array ==== */

    l_sqr  = l_size * l_size;
    large  = 0;                   /* largest cluster found. Also stack top */
    mlarge = 0;                   /* the largest cluster ever encountered  */

    for (i = 0; i < l_sqr; i++)
      {
       lptr[i] = 0;
       nsc[i]  = 0;
       vec1[i] = 0;
    }

    maxi  = l_sqr + 2;

   /* -------------------------------------------------------------------- */
   /*                                                                      */
   /*             M O N T E   C A R L O   P A R T                          */
   /*                                                                      */
   /* -------------------------------------------------------------------- */

   for (mcs = 1; mcs <= mcsmax; mcs++)
     {

      /* ==== Go through the system regularly and generate one  ==== */
      /* ==== configuartion                                     ==== */

      for (iy = 0; iy < l_size; iy++)
      {
         hold = r250( l_size,ran_ptr,mf_ptr,hold);

         for (ix = 0; ix < l_size; ix++)
         {
            if (ran[ix] < pB)
            {
               iss[ix][iy] = 1;
            } else {
               iss[ix][iy] = 0;
            }
         }
      }


         /* -------------------------------------------------------------- */
         /*                                                                */
         /*             D r o p l e t   A n a l y s i s                    */
         /*                                                                */
         /*    Free Bounary Conditions                                     */
         /*                                                                */
         /*                                                                */
         /* -------------------------------------------------------------- */

         /* ==== The array <<<row>>> holds always the previously ==== */
         /* ==== analysed row from the array iss.                ==== */
         /* ==== First row is taken as entirely unoccupied, i.e. ==== */
         /* ==== free boundary conditions                        ==== */

         cl = 0;               /* Will be the largest cluster pointer */

         for (irow = 0; irow <= l_size; irow++)
           {
            row[irow] = maxi;
         }

         /* ==== Now we go through the lattice row by row ==== */

         for (irow = 0; irow < l_size; irow++)
           {

            for (icol = 0; icol < l_size; icol++)
              {

               if (iss[irow][icol] != 1)
                 {
                  row[icol+1] = maxi;
               }
               else
		 {

		  /*  ==== see if the site is connected ==== */

                  up   = row[icol+1];
                  left = row[icol];

                  if (up != maxi)
                    {
                     /* ==== Site is connected to the previous row. ==== */

                     if (lptr[up] < 0)
                       {
                        /* ==== We found a negative label, signaling a ==== */
                        /* ==== a pointer. Search now for the root.    ==== */

                        ms = lptr[up];
                        while ( ms < 0 )
                          {
                           la = -ms;
                           ms = lptr[la];
                        }
                        lptr[up] = -la;
                        up       = la;
                     }
                  }

                  if (left != maxi)
                    {
                     /* ==== Site is connected to left neighbour. ==== */

                     if (lptr[left] < 0)
                       {
                        /* ==== We found a negative label, signaling a ==== */
                        /* ==== a pointer. Search now for the root.    ==== */

                        ms = lptr[left];
                        while ( ms < 0 )
                          {
                           la = -ms;
                           ms = lptr[la];
                        }
                        lptr[left] = -la;
                        left       = la;
                     }
                  }

                  mini = (up < left) ? up : left;
                  if (mini == maxi)
                    {
                     /* ==== Site is not connected. Assign new label ==== */

                     cl++;
                     row[icol+1] = cl;
                     iss[irow][icol] = cl;
	             lptr[cl]    = 1;

		          }
                  else
                    {
                     /* ==== Site is connected. Find minimum label ==== */

                     nofs = 1;
                     if (up != left)
                       {
                        /* ==== Possibly two clusters joined into one.  ==== */
                        /* ==== Count both parts towards the number of  ==== */
                        /* ==== in the cluster.                         ==== */

                        if (up != maxi)
                          {
                           nofs    += lptr[up];
                           lptr[up] = -mini;
                        }

                        if (left != maxi)
                          {
                           nofs      += lptr[left];
                           lptr[left] = -mini;
                        }
                     }
                     else
                       {
                         /* ==== The two possible cluster are part of  ==== */
                         /* ==== one. Count only one                   ==== */

                         nofs      += lptr[left];
                         lptr[left] = -mini;
                     }
                     row[icol+1] = mini;
                     iss[irow][icol] = mini;
                     lptr[mini]  = nofs;
                  }
               }
            }
         }

         /* -------------------------------------------------------------- */
         /*                                                                */
         /*          E n d    D r o p l e t   A n a l y s i s              */
         /*                                                                */
         /* -------------------------------------------------------------- */

         /* -------------------------------------------------------------- */
         /*                                                                */
         /*                 Analysis of droplet numbers                    */
         /*                                                                */
         /* -------------------------------------------------------------- */


         for (i = 0; i <= large; i++)  /* reset the array */
           {
            vec1[i] = (float)0.;
         }

         /* ==== Get the largest droplet and the size distribution ==== */

         large = 0;
   
	     for (i = 1; i <= cl; i++)
           {
            nis = lptr[i];
            if ( nis > 0 )
              {
               if ( nis > large ) { large = nis; largeLabel = i;}
               vec1[nis] += 1;
            }
         }

         /* ==== Accumulate ==== */


         if (large > mlarge) { mlarge = large; }

         for (i = 0; i <= large; i++)
           {
            nsc[i] += vec1[i];
         }

         /* === plot the largest cluster in the system === */

         for (iy = 0; iy < l_size; iy++)
         {

            for (ix = 0; ix < l_size; ix++)
            {
              if (iss[ix][iy] != maxi)
              {
                 if (lptr[iss[ix][iy]] < 0)
                 {
                    ms = lptr[iss[ix][iy]];
                    while ( ms < 0 )
                    {
                       la = -ms;
                       ms = lptr[la];
                    }
                    lptr[up] = -la;
                    iss[ix][iy]  = la;
                  }
              }


            }
        }
 
         /* -------------------------------------------------------------- */
         /*                                                                */
         /*               End analysis of droplet numbers                  */
         /*                                                                */
         /* -------------------------------------------------------------- */

 

   }  /* end of Monte Carlo loop */
   return 0;
} /* end of main program */
