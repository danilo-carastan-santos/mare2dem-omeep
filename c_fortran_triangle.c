// !-----------------------------------------------------------------------
// !
// !    Copyright 2009
// !    Kerry Key, David Myer
// !    Scripps Institution of Oceanography
// !    kkey@ucsd.edu
// !
// !    This file is part of MARE2DEM.
// !
// !    MARE2DEM is free software: you can redistribute it and/or modify
// !    it under the terms of the GNU General Public License as published by
// !    the Free Software Foundation, either version 3 of the License, or
// !    (at your option) any later version.
// !
// !    MARE2DEM is distributed in the hope that it will be useful,
// !    but WITHOUT ANY WARRANTY; without even the implied warranty of
// !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// !    GNU General Public License for more details.
// !
// !    You should have received a copy of the GNU General Public License
// !    along with MARE2DEM.  If not, see <http://www.gnu.org/licenses/>.
// !
// !-----------------------------------------------------------------------
//
// An interface routine to call Triangle.c from Fortran.  This aint pretty, but it seems to work.
// 
// Kerry Key and David Myer
// Scripps Institution of Oceanography
//
// Versions:
// 1.0 	March/April 2009
//
//

/* system headers */
#include <stdio.h>
#include <stdlib.h>

/* Stuff */
#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

/* local header files */
#include "triangle.h"


struct triangulateio out; 


// Helper function to clear a triangulation structure
void clear_triangulate( struct triangulateio *st )
{
    st->pointlist                = (REAL *)NULL;
    st->pointattributelist       = (REAL *)NULL;
    st->pointmarkerlist          = (int *) NULL;
    st->numberofpoints           = 0;
    st->numberofpointattributes  = 0;
    
    st->trianglelist             = (int *) NULL;
    st->triangleattributelist    = (REAL *)NULL;
    st->trianglearealist         = (REAL *)NULL;
    st->neighborlist             = (int *) NULL;
    st->numberoftriangles        = 0;
    st->numberofcorners          = 0;
    st->numberoftriangleattributes = 0;

    st->segmentlist              = (int *) NULL;
    st->segmentmarkerlist        = (int *) NULL;
    st->numberofsegments         = 0;

    st->holelist                 = (REAL *)NULL;
    st->numberofholes            = 0;

    st->regionlist               = (REAL *)NULL;
    st->numberofregions          = 0;

    st->edgelist                 = (int *) NULL;
    st->edgemarkerlist           = (int *) NULL;
    st->normlist                 = (REAL *)NULL;
    st->numberofedges            = 0;

}   // clear_triangulate

// Helper function to deallocate members of a triangulation structure
// that have been allocated.  Counts on clear_triangulate() having been
// called previous to begin of allocation so that pointers are NULL when
// not used.
void free_triangulate( struct triangulateio *st )
{
    if( st->pointlist )              free( st->pointlist );
    if( st->pointattributelist )     free( st->pointattributelist );
    if( st->pointmarkerlist )        free( st->pointmarkerlist );
    
    if( st->trianglelist )           free( st->trianglelist );
    if( st->triangleattributelist )  free( st->triangleattributelist );
    if( st->trianglearealist )       free( st->trianglearealist );
    if( st->neighborlist )           free( st->neighborlist );

    if( st->segmentlist )            free( st->segmentlist );
    if( st->segmentmarkerlist )      free( st->segmentmarkerlist );

//     if( st->holelist )               free( st->holelist );  
// 
//     if( st->regionlist )             free( st->regionlist );
    
    if( st->edgelist )               free( st->edgelist );
    if( st->edgemarkerlist )         free( st->edgemarkerlist );
    if( st->normlist )               free( st->normlist );

    // Now clear the structure so that free'd pointers are NULL
    clear_triangulate( st );

}   // free_triangulate


//--------------------------------------------------------//	
// Step 1:  Receive the grid to refine and refine with 
//          a call to Triangle.c
//
// Call this function from Fortran. Pass in the grid 
// to be refined or the starting boundary model.
// Only the refined grid array sizes are passed back to Fortran,
// where you will need to allocate new arrays 
// and then call step 2 here to pass
// the new grid arrays back. 
//--------------------------------------------------------//	


 void c_fortran_triangluate_input_(  
         char *tristr,
		 REAL *px, REAL *py , REAL *pointattributelist, int *pointmarkerlist,
		 int *numberofpoints,  int *numberofpointattributes, 
		 int *trianglelist, REAL *triangleattributelist, REAL *trianglearealist,
		 int *numberoftriangles,  int *numberofcorners, int *numberoftriangleattributes,
		 int *segmentlist, int *segmentmarkerlist, int *numberofsegments,
		 REAL *holelist,  int *numberofholes,
		 REAL *regionlist, int *numberofregions, 
		 int *onumberofpoints,  int *onumberofpointattributes,
		 int *onumberoftriangles,  int *onumberofcorners, int *onumberoftriangleattributes,
         int *onumberofsegments   ) 
{
    struct triangulateio in;
    int i;

    // Bundle up input arrays into triangulateio structure:
    clear_triangulate( &in );

    // points:
    in.numberofpoints = *numberofpoints;
    in.numberofpointattributes = *numberofpointattributes;

    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    for (i = 0; i<in.numberofpoints ; i++)  
    {
        in.pointlist[2*i] = px[i];
        in.pointlist[2*i+1] = py[i];
    }

    in.pointattributelist = (REAL *) malloc(in.numberofpoints *
                                            in.numberofpointattributes *
                                            sizeof(REAL));
    for (i = 0; i<in.numberofpoints *in.numberofpointattributes; i++)
        in.pointattributelist[i] = pointattributelist[i];

    in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
    for (i = 0; i < in.numberofpoints; i++)
        in.pointmarkerlist[i] = pointmarkerlist[i];

    // triangles:
    // if they exist then copy them over here, else just skip this section:

    in.numberoftriangles = *numberoftriangles;
    in.numberofcorners = *numberofcorners;
    in.numberoftriangleattributes = *numberoftriangleattributes;

    if (in.numberoftriangles > 0) 
    {
        // Allocate arrays
        in.trianglelist = (int *) malloc(in.numberoftriangles*in.numberofcorners* sizeof(int));
        in.triangleattributelist = (REAL *) malloc(in.numberoftriangleattributes * in.numberoftriangles  * sizeof(REAL));
        in.trianglearealist = (REAL *) malloc(in.numberoftriangles * sizeof(REAL));

        /// copy across data:	
        for (i = 0; i < in.numberoftriangles*in.numberofcorners; i++) 
            in.trianglelist[i] = trianglelist[i];
        for (i = 0; i < in.numberoftriangleattributes * in.numberoftriangles; i++) 
            in.triangleattributelist[i] = triangleattributelist[i];
        for (i = 0; i < in.numberoftriangles; i++) 
            in.trianglearealist[i] = trianglearealist[i];
    }


    // segments:	  
    in.numberofsegments = *numberofsegments;
    if (in.numberofsegments > 0) 
    {
        in.segmentlist = (int *) malloc(in.numberofsegments * 2* sizeof(int));
        in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
        for (i = 0; i < in.numberofsegments*2; i++) 
            in.segmentlist[i] = segmentlist[i];
        for (i = 0; i < in.numberofsegments; i++) 
            in.segmentmarkerlist[i] = segmentmarkerlist[i];
    }

    // holes:	  
    in.numberofholes = *numberofholes;
    if (in.numberofholes > 0) 
    {
        in.holelist = (REAL *) malloc(in.numberofholes * 2 * sizeof(REAL));
        for (i = 0; i < in.numberofholes*2; i++) 
            in.holelist[i] = holelist[i];
    }
    // regions:	  
    in.numberofregions = *numberofregions;
    if (in.numberofregions > 0) 
    {
        in.regionlist = (REAL *) malloc(in.numberofregions * 4 * sizeof(REAL));
        for (i = 0; i < in.numberofregions*4; i++) 
            in.regionlist[i] = regionlist[i];
    }

    // edges:	 
    // 	  in.numberofedges = *numberofedges;
    // 	  in.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    // 	  in.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
    // 	  in.normlist = (REAL *) NULL;;
    // 	  



    /* Make necessary initializations so that Triangle can return a */
    /*   triangulation in `out' */
    clear_triangulate( &out );


    // Triangulate that mofo:
    triangulate(tristr, &in, &out, (struct triangulateio *) NULL);


    // copy array sizes to send back to Fortran:
    *onumberofpoints = out.numberofpoints;
    *onumberofpointattributes  = out.numberofpointattributes ;
    *onumberoftriangles = out.numberoftriangles;
    *onumberofcorners = out.numberofcorners;
    *onumberoftriangleattributes = out.numberoftriangleattributes;
    *onumberofsegments = out.numberofsegments;

    // New grid dimensions are then returned to FORTRAN

    // Free up allocated memory	for input mesh:
    
   if ( in.holelist )   free( in.holelist );  //kwk nov 2011: moved this before free_triangulate to avoid memory leak, aack!
   if ( in.regionlist ) free( in.regionlist );
   
    free_triangulate( &in );

}

//--------------------------------------------------------//	
// Step 2:  Pass back the refined grid to Fortran:
//
// Function for outputing the new grid back to Fortran.
// after Fortran has allocated space for the new arrays
//--------------------------------------------------------//
void c_fortran_triangluate_output_( 
   REAL *px, REAL *py, REAL *pointattributelist, int *pointmarkerlist,
   int *trianglelist, REAL *triangleattributelist,  
   int *neighborlist,
   int *segmentlist, int *segmentmarkerlist ) 
{
    int i;

    // Copy over the arrays:
    // points:
    //  printf("copying points \n");
    for (i = 0; i < out.numberofpoints; i++) { px[i] = out.pointlist[2*i]; py[i] = out.pointlist[2*i+1]; } ;
    for (i = 0; i < out.numberofpoints*out.numberofpointattributes; i++)	pointattributelist[i] = out.pointattributelist[i];
    for (i = 0; i < out.numberofpoints; i++) pointmarkerlist[i] = out.pointmarkerlist[i];

    // triangles:
    // printf("copying triangles, #tris = %i\n",out.numberoftriangles);
    if (out.numberoftriangles > 0 ) {
        for (i = 0; i < out.numberoftriangles*out.numberofcorners; i++)  trianglelist[i] = out.trianglelist[i];
        for (i = 0; i < out.numberoftriangles * 3; i++) neighborlist[i] = out.neighborlist[i];
        if (out.numberoftriangles*out.numberoftriangleattributes > 0 ) {
            for (i = 0; i < out.numberoftriangleattributes * out.numberoftriangles; i++) 
                triangleattributelist[i] = out.triangleattributelist[i] ;
        }	   
    }



    // segments:
    if (out.numberofsegments > 0 ) {
        //  printf("copying segments \n");
        for (i = 0; i < out.numberofsegments*2; i++)  segmentlist[i] = out.segmentlist[i];
        for (i = 0; i < out.numberofsegments; i++) segmentmarkerlist[i] = out.segmentmarkerlist[i];
    }

    // edges:	 
    // 	  if (out.numberofedges > 0 ) {
    // 	  printf("copying edges \n");
    // 	  for (i = 0; i < out.numberofedges*2; i++)  edgelist[i] = out.edgelist[i];             /* Needed only if -e switch used. */
    // 	  for (i = 0; i < out.numberofedges; i++) edgemarkerlist[i] = out.edgemarkerlist[i];   /* Needed if -e used and -B not used. */
    // 	   printf("done with edges , = \n");
    // 	  }


    // Free up allocated memory for output mesh:
    free_triangulate( &out );
   
}
