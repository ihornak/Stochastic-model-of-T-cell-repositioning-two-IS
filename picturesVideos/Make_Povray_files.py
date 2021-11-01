#!/usr/bin/env python2
import os,sys,glob,copy,re,shutil,collections,string
import glob
import re
from os.path import join, dirname



def cytoskeleton( name_of_the_file , MTOC_file , IS_catching_name , IS_capture_shrinkage , IS_capture_shrinkage_2 , IS_cortical_slid , IS_cortical_slid_2 ,
 Dynein_surface_random_dist ,
 dynein_abscissa , dynein_abscissa_attachment , Dynein_IS_capt , Dynein_IS_capt_2 , nucleus_param ):



    # loading of all microtubulas
    file = open( name_of_the_file , 'r')
    microtubule_TMP = file.read()
    list_of_microtubules = microtubule_TMP.split('END')
    file.close()
    if( list_of_microtubules[ - 1 ] == "\n" ):
        list_of_microtubules.pop()




    file = open( MTOC_file , 'r')
    MTOC_TMP = file.read()
    list_MTOC_polygon = MTOC_TMP.split('END')
    file.close()
    list_MTOC_polygon.pop();


    file = open( IS_catching_name , 'r')
    IS_catching_TMP = file.read()
    list_IS_catching_TMP = IS_catching_TMP.split('END')
    file.close()
    list_IS_catching_TMP.pop();

 







    file = open( Dynein_surface_random_dist , 'r')
    Dynein_surface_random_TMP = file.read()
    Dynein_surface_random_list = Dynein_surface_random_TMP.split('\nEND\n')
    file.close()
    Dynein_surface_random_list.pop()

    file = open( Dynein_IS_capt , 'r')
    Dynein_IS_capt_TMP = file.read()
    Dynein_IS_capt_list = Dynein_IS_capt_TMP.split('\nEND\n')
    file.close()
    Dynein_IS_capt_list.pop()


    file_2 = open( Dynein_IS_capt_2 , 'r')
    Dynein_IS_capt_TMP_2 = file_2.read()
    Dynein_IS_capt_list_2 = Dynein_IS_capt_TMP_2.split('\nEND\n')
    file_2.close()
    Dynein_IS_capt_list_2.pop()





    file = open( dynein_abscissa , 'r')
    dynein_abscissa_TMP = file.read()
    dynein_abscissa_list = dynein_abscissa_TMP.split('\nEND\n')
    file.close()
    dynein_abscissa_list.pop();


    file = open( dynein_abscissa_attachment , 'r')
    dynein_abscissa_attachment_TMP = file.read()
    dynein_abscissa_attachment_list = dynein_abscissa_attachment_TMP.split('\nEND\n')
    file.close()
    dynein_abscissa_attachment_list.pop();





    #adjusting points of MTOC to the form of one list with elements composed of vectors
    list_of_MTOC_points = []
    for polygon in list_MTOC_polygon:
        tmp = polygon.split('\n')
        #dispose '' at the end
        result = list( filter( None , tmp ) )
        list_of_MTOC_points.extend( result )


    MTOC_center = list_of_MTOC_points[ 0 ]


    MTOC_SPHERE = """sphere {
<%s >, 0.4
texture{ pigment{color Black transmit 0.1}   finish {phong 0.2}}
}
""" % (MTOC_center )


    list_IS_catching = []
    for polygon in list_IS_catching_TMP:
        tmp = polygon.split('\n')
        #dispose '' at the end
        result = list( filter( None , tmp ) )
        list_IS_catching.extend( result )
    # all microtubulas are loaded and saved into list_of_microtubules


    # string adjusting camera, background and light source
    # it will be included just once   
    includes = """
    #include \"colors.inc\"
    #include \"strings.inc\"
    #include "shapes.inc"
    #include "stones.inc"

    background
    {
    color rgb <1,1,1>
    }
    camera
    {
        
        location < 0 , 12.0 , 0 >   sky   <0,0,1> look_at < 0.0 , 0.0 , 0.0 >  
        }
        light_source { < - 4, 4 ,- 2 > 1 shadowless }
        """

    #I load here parametres of the cell
    file_cell = open( 'cell_shapes/Cell_Parametres.txt' , 'r')
    numbers_cell = file_cell.read()
    axis_2 = numbers_cell.split(',')
    file_cell.close()


    ellipsoid = """
    object{
    Spheroid( //CenterVector,
             <0.0,0.0,0.00>,
             // RadiusVector Rx,Ry,Rz )
             <%s,%s,%s> )
    texture{ pigment{color rgbt<.75,.2,0,.7>  transmit 0.8}
            finish { phong 1}
          } // end of texture
    //scale<1,1,1>
    rotate<0,0,0>
    translate<0,0.0,0>
    no_shadow
    }"""%( axis_2[ 0 ] , axis_2[ 0 ] , axis_2[ 1 ]  )



    IS_lower_part_of_cell = """cylinder{ <0.0 , 0.0 , 0.0 >,<0.0 , 0.0 , -5.3 >,4.5
          texture{ pigment{ color rgb <0.9921569,  0.5686275,  0.2823529> transmit 0.8 }
                   finish { phong 1}
                 }
        }"""


    #I load here parametres of the cells nucleus
    file_cell_nucleus = open( 'cell_shapes/Nucleus_Parametres.txt' , 'r')
    numbers_cell_nucleus = file_cell_nucleus.read()
    axis_and_center = numbers_cell_nucleus.split('\nEND\n')


    axis_2_nucleus = axis_and_center[ 0 ].split(',')
    center_nucleus = axis_and_center[ 1 ].split(',')
    file_cell_nucleus.close()



    file = open( nucleus_param , 'r')
    nucleus_file = file.read()
    nucleus_data = nucleus_file.split('\n')
    file.close()
    nucleus_data.pop()
    nucleus_data.pop()
    center_nucleus = nucleus_data[ 0 ].split(',')
    axis = nucleus_data[ 1 ].split(',')

    nucleus_ellipsoid = """
    object{
    Spheroid( //CenterVector,
             <%s,%s,%s>,
             // RadiusVector Rx,Ry,Rz )
             <%s,%s,%s> )
    texture{ pigment{color rgb <  0.2470588,  0.4431373,  0.9058824  > transmit 0.5}
            finish { phong 1}
          } // end of texture
    //scale<1,1,1>
    rotate<0,0,0>
    translate<0,0.0,0>
    no_shadow
    }"""%( center_nucleus[ 0 ] , center_nucleus[ 1 ] , center_nucleus[ 2 ] , axis[ 0 ] , axis[ 1 ] , axis[ 1 ] )







    sphere = """
    sphere {
    <0,0,0.0>, 3.0
    texture{ pigment{color Gray transmit 0.9}   finish {phong 0.2}}
    }
    """

    nucleus = """
    sphere {
    <0,0,0.0>, 1.6
    texture{ pigment{color Gray transmit 0.9}
    finish {phong 0.2}}
    }
    """



    list_of_splines = []
    for microtubule in range(0,len( list_of_microtubules )):        # Second Example
        spline_template = """
        #declare MicrotubuleSpline_%d =
        spline {
        linear_spline %s 
        }
        """% (
        microtubule, list_of_microtubules[ microtubule ])
        list_of_splines.append( spline_template )



    dynein_index = []
    for microtubule in range(0,len( list_of_microtubules )):
        kkk = list_of_microtubules[ microtubule ].split('>')
        aaa = kkk[ -1 ]
        aaa = aaa[1:]
        aaa = aaa[:-1]
        list_of_microtubules[ microtubule ] = list_of_microtubules[ microtubule ].rstrip()
        dynein_index.append(aaa)
        list_of_microtubules[ microtubule ] = list_of_microtubules[ microtubule ][:-2]



    list_of_splines = []
    for microtubule in range(0,len( list_of_microtubules )):        # Second Example
        spline_template = """
        #declare MicrotubuleSpline_%d =
        spline {
        cubic_spline %s
        }
        """% (
        microtubule, list_of_microtubules[ microtubule ])
        list_of_splines.append( spline_template )


    list_of_print_strings = [ ]
    for microtubule in range( 0,len( list_of_microtubules )):        # Second Example
        if int( dynein_index[ microtubule ] ) == 0:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.02
        pigment { Blue }
        no_shadow
        }
        #declare ctr = ctr + 0.0001;
        #end

        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 1:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.02
        pigment { Blue }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 2:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.02
        pigment { Green }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 3:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.02
        pigment { Black }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 4:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.015
        pigment { Black }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 9:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.03
        pigment { Yellow }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 10:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.03
        pigment { Yellow }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 11:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.03
        pigment { Yellow }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 20:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.03
        pigment { Red }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )
        if int( dynein_index[ microtubule ] ) == 40:
            print_string = """
        #declare ctr = 0;
        #while (ctr < 2)
        sphere {
        MicrotubuleSpline_%d(ctr),.03
        pigment { Black }
        }
        #declare ctr = ctr + 0.0001;
        #end
        """% (
        microtubule)
            list_of_print_strings.append( print_string )







    list_of_MTOC_SPHERES = []
    d = 0
    for MTOC_point in list_of_MTOC_points:        # Second Example
        a = ( d - 1 ) / 5
        if d== 0:
            color_string = "Black"
        if a == 0:
            color_string = "Blue"
        elif a == 1:
            color_string = "Red"
        elif a == 2:
            color_string = "Green"
        elif a == 3:
            color_string = "Aquamarine"
        spline_template = """sphere {
<%s >, 0.04
texture{ pigment{color %s transmit 0.1}   finish {phong 0.2}}
}
""" % (MTOC_point , color_string )
        list_of_MTOC_SPHERES.append( spline_template)
        d = d + 1



    list_of_IS_catching_SPHERES = []
    for IS_catching in list_IS_catching:        # Second Example
        spline_template = """sphere {
<%s >, 0.1
texture{ pigment{color Black transmit 0.0}   finish {phong 0.2}}
}
""" % (IS_catching )# , color_string
        list_of_IS_catching_SPHERES.append( spline_template)










    Dynein_surface_random_list_SPHERES = []
    for point in Dynein_surface_random_list:        # Second Example
                spline_template = """sphere {
        <%s >, 0.017
        texture{ pigment{color LimeGreen transmit 0.1}   finish {phong 0.2}}
        }
        """ % ( point )# , color_string
                Dynein_surface_random_list_SPHERES.append( spline_template )


    Dynein_IS_capt_list_SHPHERE = []
    for point in Dynein_IS_capt_list:        # Second Example
                spline_template = """sphere {
        <%s >, 0.05
        texture{ pigment{color BlueViolet transmit 0.1}   finish {phong 0.2}}
        }
        """ % ( point )# , color_string
                Dynein_IS_capt_list_SHPHERE.append( spline_template )



    Dynein_IS_capt_list_SHPHERE_2 = []
    for point in Dynein_IS_capt_list_2:        # Second Example
                spline_template = """sphere {
        <%s >, 0.05
        texture{ pigment{color BlueViolet transmit 0.1}   finish {phong 0.2}}
        }
        """ % ( point )# , color_string
                Dynein_IS_capt_list_SHPHERE_2.append( spline_template )






    dynein_abscissa_list_SPHERES = []
    for point in dynein_abscissa_list:        # Second Example
                spline_template = """sphere {
        <%s >, 0.11
        texture{ pigment{color Black transmit 0.1}   finish {phong 0.2}}
        }
        """ % ( point )# , color_string
                dynein_abscissa_list_SPHERES.append( spline_template )


    dynein_abscissa_attachment_list_SPHERES = []
    for point in dynein_abscissa_attachment_list:        # Second Example
                spline_template = """sphere {
        <%s >, 0.1
        texture{ pigment{color VeryDarkBrown transmit 0.1}   finish {phong 0.2}}
        }
        """ % ( point )# , color_string
                dynein_abscissa_attachment_list_SPHERES.append( spline_template )







    file_capture_shrinkage = open( IS_capture_shrinkage , 'r')
    capture_shrinkage_content = file_capture_shrinkage.read()
    capture_shrinkage_content = capture_shrinkage_content[ : -2 ]
    capture_shrinkage_content_2 = capture_shrinkage_content.split('\nEND\n')
    file_capture_shrinkage.close()

    IS_capture_shrinkage = """cylinder{ <%s>,<%s>,%s
          texture{ pigment{ color Coral  transmit 0.52}
                   finish { phong 1}
                 }
        }"""% (capture_shrinkage_content_2[ 0 ] , capture_shrinkage_content_2[ 1 ] ,  capture_shrinkage_content_2[ 2 ]  )


    file_capture_shrinkage_b = open( IS_capture_shrinkage_2 , 'r')
    capture_shrinkage_content_b = file_capture_shrinkage_b.read()
    capture_shrinkage_content_b = capture_shrinkage_content_b[ : -2 ]
    capture_shrinkage_content_2_b = capture_shrinkage_content_b.split('\nEND\n')
    file_capture_shrinkage_b.close()

    IS_capture_shrinkage_2 = """cylinder{ <%s>,<%s>,%s
          texture{ pigment{ color Coral  transmit 0.52}
                   finish { phong 1}
                 }
        }"""% (capture_shrinkage_content_2_b[ 0 ] , capture_shrinkage_content_2_b[ 1 ] ,  capture_shrinkage_content_2_b[ 2 ]  )





    file_IS_cortical_slid = open( IS_cortical_slid , 'r')
    IS_cortical_slid_content = file_IS_cortical_slid.read()
    IS_cortical_slid_content = IS_cortical_slid_content[ : -2 ]
    IS_cortical_slid_content_2 = IS_cortical_slid_content.split('\nEND\n')
    file_IS_cortical_slid.close()

    IS_cortical_slid = """cylinder{ <%s>,<%s>,%s
          texture{ pigment{ color Cyan   transmit 0.5}
                   finish { phong 1}
                 }
        }"""% (IS_cortical_slid_content_2[ 0 ] , IS_cortical_slid_content_2[ 1 ] ,  IS_cortical_slid_content_2[ 2 ]  )



    file_IS_cortical_slid_b = open( IS_cortical_slid_2 , 'r')
    IS_cortical_slid_content_b = file_IS_cortical_slid_b.read()
    IS_cortical_slid_content_b = IS_cortical_slid_content_b[ : -2 ]
    IS_cortical_slid_content_2_b = IS_cortical_slid_content_b.split('\nEND\n')
    file_IS_cortical_slid_b.close()

    IS_cortical_slid_b = """cylinder{ <%s>,<%s>,%s
          texture{ pigment{ color Cyan   transmit 0.5}
                   finish { phong 1}
                 }
        }"""% (IS_cortical_slid_content_2_b[ 0 ] , IS_cortical_slid_content_2_b[ 1 ] ,  IS_cortical_slid_content_2_b[ 2 ]  )




    # creating name of the file
    #taking time from name
    b = re.split('[_ .]',name_of_the_file)
    name_of_file = 'microtubule_'+b[2] + '.pov'
    # opening file
    save_file = "./povrayFiles/" + name_of_file
    f = open( save_file , 'w')
    f.write( includes )

    f.write( ellipsoid )
    f.write( nucleus_ellipsoid )


    f.write( IS_capture_shrinkage )
    f.write( IS_cortical_slid )

    f.write( IS_capture_shrinkage_2 )
    f.write( IS_cortical_slid_b )
    f.write( MTOC_SPHERE )
    #f.write( IS_lower_part_of_cell )

    
    for microtubule in range(0,len( list_of_splines )):
        f.write( list_of_splines[ microtubule ] )

    for microtubule in range(0,len( list_of_print_strings )):
        f.write( list_of_print_strings[ microtubule ] )

    #for MTOC_point in list_of_MTOC_SPHERES:
        #f.write( MTOC_point )

    #capture-shrinkage dyneins
    for IS_catching in list_of_IS_catching_SPHERES:
        f.write( IS_catching )

    #for sphere in Dynein_surface_random_list_SPHERES:
                #f.write( sphere )

    #for sphere in Dynein_IS_capt_list_SHPHERE:
                #f.write( sphere )


    #for sphere in Dynein_IS_capt_list_SHPHERE_2:
                #f.write( sphere )


    for sphere in dynein_abscissa_list_SPHERES:
                f.write( sphere )
                
    #for sphere in dynein_abscissa_attachment_list_SPHERES:
                #f.write( sphere )



    f.close()
    return







list_of_files = glob.glob("./textFiles/microtubule*")
number_of_files_tmp = len( list_of_files )

    

for number  in range( 0 , number_of_files_tmp ):
    micro_file = """./textFiles/microtubule_%s.txt"""%(
        number)
    MTOC_file = """./textFiles/MTOC_%s.txt"""%(
        number)
    IS_catching_file = """./textFiles/IS_catching_%s.txt"""%(
        number )

    IS_cortical_sl = """./textFiles/IS_cortical_sl_%s.txt"""%(
        number )

    IS_capture_shrinkage = """./textFiles/IS_capture_shrinkage_0.txt"""


    Dynein_surface_random_dist = """./textFiles/Dynein_surface_randomly_distributed_%s.txt"""%(
        number )
    Dynein_IS_capt = """./textFiles/Dynein_IS_capture_%s.txt"""%(
        number )
    dynein_abscissa = """./textFiles/dynein_abscissa_%s.txt"""%(
            number )

    dynein_abscissa_attachment = """./textFiles/dynein_abscissa_attachment_%s.txt"""%(
            number )


    IS_cortical_sl_2 = """./textFiles/IS_cortical_sl_2_0.txt"""
    

    IS_capture_shrinkage_2 = """./textFiles/IS_capture_shrinkage_2_%s.txt"""%(
        number )

    Dynein_IS_capt_2 = """./textFiles/Dynein_IS_capture_2_%s.txt"""%(
        number )


    nucleus_param = """./textFiles/nucleus_print_%s.txt"""%(
        number )



    cytoskeleton( micro_file , MTOC_file  , IS_catching_file , IS_capture_shrinkage , IS_capture_shrinkage_2 , IS_cortical_sl , IS_cortical_sl_2 , Dynein_surface_random_dist 
    , dynein_abscissa , dynein_abscissa_attachment , Dynein_IS_capt , Dynein_IS_capt_2 , nucleus_param )


