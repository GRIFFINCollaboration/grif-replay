
////////////////////////////
// Helper functions for GRIFFIN and QED geometry
////////////////////////////

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "grif-angles.h"

float acos_table[ACOS_TABLE_SIZE]; // Global acos lookup table

// calculate dot product of two three-dimensional vectors
double dot_product(double vector1[], double vector2[]) {
  double result = 0.0;
  result += vector1[0] * vector2[0];
  result += vector1[1] * vector2[1];
  result += vector1[2] * vector2[2];
  return result;
}

// calculate cross product of two three-dimensional vectors
void cross_product(double vector1[], double vector2[], double result[]) {
  result[0] = vector1[1] * vector2[2] - vector1[2] * vector2[1];
  result[1] = vector1[2] * vector2[0] - vector1[0] * vector2[2];
  result[2] = vector1[0] * vector2[1] - vector1[1] * vector2[0];
}

// calculate the product of the magnitudes of two three-dimensional vector
double vector_magnitude_product(double vector1[], double vector2[]){
  double result1 = 0.0,result2 = 0.0;
  result1 += vector1[0] * vector1[0];
  result1 += vector1[1] * vector1[1];
  result1 += vector1[2] * vector1[2];
  result1 = sqrt(result1);
  result2 += vector2[0] * vector2[0];
  result2 += vector2[1] * vector2[1];
  result2 += vector2[2] * vector2[2];
  result2 = sqrt(result2);

  return (result1*result2);
}

// rotate a three-dimensional vector, v, around an axis of unit vector (u) by the given angle
void rotate_vector(double v[], double u[], double angle, double result[]) {
  double dot,cross[3];
  dot = dot_product(u,v);
  cross_product(u,v,cross);
  result[0] = v[0]*cos(DEGREES_TO_RADIANS*angle) + cross[0]*sin(DEGREES_TO_RADIANS*angle) + u[0]*dot*(1-cos(DEGREES_TO_RADIANS*angle));
  result[1] = v[1]*cos(DEGREES_TO_RADIANS*angle) + cross[1]*sin(DEGREES_TO_RADIANS*angle) + u[1]*dot*(1-cos(DEGREES_TO_RADIANS*angle));
  result[2] = v[2]*cos(DEGREES_TO_RADIANS*angle) + cross[2]*sin(DEGREES_TO_RADIANS*angle) + u[2]*dot*(1-cos(DEGREES_TO_RADIANS*angle));
  if(result[0]>-0.0001 && result[0]<0.0001){ result[0]=0.0; }
  if(result[1]>-0.0001 && result[1]<0.0001){ result[1]=0.0; }
  if(result[2]>-0.0001 && result[2]<0.0001){ result[2]=0.0; }
}

// Initialize table mapping index to acos(x)
void init_acos_table(void){
  for(int i = 0; i < ACOS_TABLE_SIZE; i++){
    // Map index [0, 511] to input range [-1.0, 1.0]
    float x = -1.0f + (2.0f * i) / (ACOS_TABLE_SIZE - 1);
    acos_table[i] = acosf(x);
  }
}

// Fast lookup function with linear interpolation because acos is slow
float fast_acos(float x){
  // Bound input to valid acos range [-1.0, 1.0]
  if(x < -1.0f) x = -1.0f;
  if(x > 1.0f) x = 1.0f;

  // Map x to table index scale [0.0, 511.0]
  float val = (x + 1.0f) * 0.5f * (ACOS_TABLE_SIZE - 1);

  int index = (int)val;
  if(index >= ACOS_TABLE_SIZE - 1){
    return acos_table[ACOS_TABLE_SIZE - 1];
  }

  // Linear interpolation for higher accuracy
  float fraction = val - index;
  float y0 = acos_table[index];
  float y1 = acos_table[index + 1];

  return y0 + fraction * (y1 - y0);
}

// Calculate the scattering angle from a HPGe energy (keV) assuming it is
// the secondary (final) photon from a single Compton scatter of an initial_energy gamma ray
int compton_angle(float ecal, float initial_energy){
  //  return (int)( RADIANS_TO_DEGREES*fast_acos( 1 - (511.0/ecal) + (511.0/initial_energy) ) );
  float factor;
  factor = 1 - (511.0/ecal) + (511.0/initial_energy);
  if(factor<-1.0 || factor>1.0){ return 0; }
  return (int)( RADIANS_TO_DEGREES*fast_acos( factor ) );
}

// Calculate the secondary photon energy from a scattering angle and an initial_energy gamma ray
int secondary_energy(float angle, float initial_energy){
  return (int)(initial_energy/(1+((initial_energy/511.0)*(1-cos(angle*DEGREES_TO_RADIANS)))));
}

// Calculate angular difference between two HPGe
double angular_diff_GeGe(int c1, int c2, int distance){
  double vec1[3], vec2[3], dot, mag, angle;

  if(distance==110){
    vec1[0] = grif_crystal_cartesian_110mm[c1][0]; vec1[1] = grif_crystal_cartesian_110mm[c1][1]; vec1[2] = grif_crystal_cartesian_110mm[c1][2];
    vec2[0] = grif_crystal_cartesian_110mm[c2][0]; vec2[1] = grif_crystal_cartesian_110mm[c2][1]; vec2[2] = grif_crystal_cartesian_110mm[c2][2];
  }else if(distance==145){
    fprintf(stderr,"angular_diff_GeGe function error: 145mm distance case not implemented yet\n");
    return 0;
    //  vec1[0] = grif_crystal_cartesian_145mm[c1][0]; vec1[1] = grif_crystal_cartesian_145mm[c1][1]; vec1[2] = grif_crystal_cartesian_145mm[c1][2];
    //  vec2[0] = grif_crystal_cartesian_145mm[c2][0]; vec2[1] = grif_crystal_cartesian_145mm[c2][1]; vec2[2] = grif_crystal_cartesian_145mm[c2][2];
  }else{
    fprintf(stderr,"angular_diff_GeGe function error: unknown distance = %d. Expected int of 110 or 145\n",distance);
    return 0;
  }
  dot = dot_product(vec1,vec2);
  mag = vector_magnitude_product(vec1,vec2);
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Calculate angular difference between two QED pixels
// pos is DSSD number [1-6], qed is pixel number [0-1023]
double angular_diff_QEDQED(int pos1, int qed1, int pos2, int qed2){
  double vec1[3], vec2[3], dot, mag, angle;

  pos1--;  pos2--;
  vec1[0] = qed_cartesian[pos1][qed1][0]; vec1[1] = qed_cartesian[pos1][qed1][1]; vec1[2] = qed_cartesian[pos1][qed1][2];
  vec2[0] = qed_cartesian[pos2][qed2][0]; vec2[1] = qed_cartesian[pos2][qed2][1]; vec2[2] = qed_cartesian[pos2][qed2][2];

  dot = dot_product(vec1,vec2);
  mag = vector_magnitude_product(vec1,vec2);
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Calculate angular difference between two QED pixels
// pos is DSSD number [1-6], qed is pixel number [0-1023]
double angular_diff_QEDGe(int pos1, int qed1, int c2, int distance){
  double vec1[3], vec2[3], dot, mag, angle;

  pos1--;
  vec1[0] = qed_cartesian[pos1][qed1][0]; vec1[1] = qed_cartesian[pos1][qed1][1]; vec1[2] = qed_cartesian[pos1][qed1][2];
  vec2[0] = grif_crystal_cartesian_110mm[c2][0]; vec2[1] = grif_crystal_cartesian_110mm[c2][1]; vec2[2] = grif_crystal_cartesian_110mm[c2][2];

  dot = dot_product(vec1,vec2);
  mag = vector_magnitude_product(vec1,vec2);
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Calculate the azimuthal angle between a scattering event in one clover and a coincidence plane defined by two HPGe
// c1 and c2 are the same as for angular correlations. These define the coincidence plane.
// c3 is another crystal in one of the clovers.
double azimuthal_GeGeGe(int c1, int c2, int c3, int distance){
  double vec1[3], vec2[3], vec3[3], coincidence_plane[3], scattering_plane[3], dot, mag, angle;

  if(distance==110){
    vec1[0] = grif_crystal_cartesian_110mm[c1][0]; vec1[1] = grif_crystal_cartesian_110mm[c1][1]; vec1[2] = grif_crystal_cartesian_110mm[c1][2];
    vec2[0] = grif_crystal_cartesian_110mm[c2][0]; vec2[1] = grif_crystal_cartesian_110mm[c2][1]; vec2[2] = grif_crystal_cartesian_110mm[c2][2];
    vec3[0] = grif_crystal_cartesian_110mm[c3][0]; vec3[1] = grif_crystal_cartesian_110mm[c3][1]; vec3[2] = grif_crystal_cartesian_110mm[c3][2];
  }else if(distance==145){
    fprintf(stderr,"azimuthal_GeGeGe function error: 145mm distance case not implemented yet\n");
    return 0;
    //  vec1[0] = grif_crystal_cartesian_145mm[c1][0]; vec1[1] = grif_crystal_cartesian_145mm[c1][1]; vec1[2] = grif_crystal_cartesian_145mm[c1][2];
    //  vec2[0] = grif_crystal_cartesian_145mm[c2][0]; vec2[1] = grif_crystal_cartesian_145mm[c2][1]; vec2[2] = grif_crystal_cartesian_145mm[c2][2];
    //  vec3[0] = grif_crystal_cartesian_145mm[c3][0]; vec3[1] = grif_crystal_cartesian_145mm[c3][1]; vec3[2] = grif_crystal_cartesian_145mm[c3][2];
  }else{
    fprintf(stderr,"azimuthal_GeGeGe function error: unknown distance = %d. Expected int of 110 or 145\n",distance);
    return 0;
  }
  cross_product(vec1,vec2,coincidence_plane); // the Normal vector of the plane (c1,c2)
  cross_product(vec2,vec3,scattering_plane);  // the Normal vector of the plane (c2,c3)
  // Now find the angle between the coincidence and scattering planes, the azimuthal
  dot = dot_product(coincidence_plane,scattering_plane);
  mag = vector_magnitude_product(coincidence_plane,scattering_plane);
  if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Calculate Compton scattering angle between one QED pixel and one HPGe from their cartesian coordinates
// pos is DSSD number [1-6], qed is pixel number [0-1023], ge is crystal number [0-63]
double scattering_angle_QEDGe(int pos, int qed, int ge){
  double vec1[3], vec2[3], vec3[3], dot, mag, angle;

  pos--; // pos is now 0-5 within this function
  vec1[0] = qed_cartesian[pos][qed][0]; vec1[1] = qed_cartesian[pos][qed][1]; vec1[2] = qed_cartesian[pos][qed][2];
  // vec2[0] = grif_crystal_cartesian_110mm[ge][0]; vec2[1] = grif_crystal_cartesian_110mm[ge][1]; vec2[2] = grif_crystal_cartesian_110mm[ge][2];
  // vec1 is the vector from origin to DSSD pixel
  // vec2 is the vector from origin to HPGe crystal
  // The dot product of vec1 and vec2 would give the angular difference between these - as required for angular correlations
  // Here we want the Compton scattering angle so we want the dot product of vec1 and vec3 where vec3 passes through the DSSD pixel and HPGe crystal
  vec3[0] = grif_crystal_cartesian_110mm[ge][0] - qed_cartesian[pos][qed][0];
  vec3[1] = grif_crystal_cartesian_110mm[ge][1] - qed_cartesian[pos][qed][1];
  vec3[2] = grif_crystal_cartesian_110mm[ge][2] - qed_cartesian[pos][qed][2];

  dot = dot_product(vec1,vec3);
  mag = vector_magnitude_product(vec1,vec3);
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Calculate Compton scattering angle between one HPGe and one QED pixel from their cartesian coordinates
// pos is DSSD number [1-6], qed is pixel number [0-1023], ge is crystal number [0-63]
double scattering_angle_GeQED(int pos, int qed, int ge){
  double vec1[3], vec2[3], vec3[3], dot, mag, angle;

  pos--; // pos is now 0-5 within this function
  //vec1[0] = qed_cartesian[pos][qed][0]; vec1[1] = qed_cartesian[pos][qed][1]; vec1[2] = qed_cartesian[pos][qed][2];
  vec2[0] = grif_crystal_cartesian_110mm[ge][0]; vec2[1] = grif_crystal_cartesian_110mm[ge][1]; vec2[2] = grif_crystal_cartesian_110mm[ge][2];
  // vec1 is the vector from origin to DSSD pixel
  // vec2 is the vector from origin to HPGe crystal
  // The dot product of vec1 and vec2 would give the angular difference between these - as required for angular correlations
  // Here we want the Compton scattering angle so we want the dot product of vec2 and vec3 where vec3 passes through the DSSD pixel and HPGe crystal
  vec3[0] = qed_cartesian[pos][qed][0] - grif_crystal_cartesian_110mm[ge][0];
  vec3[1] = qed_cartesian[pos][qed][1] - grif_crystal_cartesian_110mm[ge][1];
  vec3[2] = qed_cartesian[pos][qed][2] - grif_crystal_cartesian_110mm[ge][2];

  dot = dot_product(vec2,vec3);
  mag = vector_magnitude_product(vec2,vec3);
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Calculate Compton scattering angle between two QED pixels from their cartesian coordinates
// pos is DSSD number [1-6], qed is pixel number [0-1023]
double scattering_angle_QEDQED(int pos1, int qed1, int pos2, int qed2){
  double vec1[3], vec2[3], vec3[3], dot, mag, angle;

  pos1--; pos2--; // pos is now 0-5 within this function
  vec1[0] = qed_cartesian[pos1][qed1][0]; vec1[1] = qed_cartesian[pos1][qed1][1]; vec1[2] = qed_cartesian[pos1][qed1][2];
  // vec2[0] = grif_crystal_cartesian_110mm[ge][0]; vec2[1] = grif_crystal_cartesian_110mm[ge][1]; vec2[2] = grif_crystal_cartesian_110mm[ge][2];
  // vec1 is the vector from origin to DSSD pixel
  // vec2 is the vector from origin to HPGe crystal
  // The dot product of vec1 and vec2 would give the angular difference between these - as required for angular correlations
  // Here we want the Compton scattering angle so we want the dot product of vec1 and vec3 where vec3 passes through the DSSD pixel and HPGe crystal
  vec3[0] = qed_cartesian[pos2][qed2][0] - qed_cartesian[pos1][qed1][0];
  vec3[1] = qed_cartesian[pos2][qed2][1] - qed_cartesian[pos1][qed1][1];
  vec3[2] = qed_cartesian[pos2][qed2][2] - qed_cartesian[pos1][qed1][2];

  dot = dot_product(vec1,vec3);
  mag = vector_magnitude_product(vec1,vec3);
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  return angle;
}

// Double Compton Scatter (DCS) - delta phi, angle between the two scattering planes
// Calculate the azimuthal angle between the two scattering planes defined by two QED-HPGe scatter events
// c1 and c2 are the QED pixel and HPGe of one event. These define the scattering plane.
// c3 and c4 are the QED pixel and HPGe of one event. These define the scattering plane.
double azimuthal_DCS(int pos1, int qed1, int ge1, int pos2, int qed2, int ge2){
  double vec1[3], vec2[3], vec3[3], vec4[3], first_scattering_plane[3], second_scattering_plane[3], dot, mag, angle;

  pos1--; // pos1 is now 0-5 within this function
  pos2--; // pos2 is now 0-5 within this function
  vec1[0] = qed_cartesian[pos1][qed1][0];         vec1[1] = qed_cartesian[pos1][qed1][1];         vec1[2] = qed_cartesian[pos1][qed1][2];
  vec2[0] = grif_crystal_cartesian_110mm[ge1][0]; vec2[1] = grif_crystal_cartesian_110mm[ge1][1]; vec2[2] = grif_crystal_cartesian_110mm[ge1][2];

  vec3[0] = qed_cartesian[pos2][qed2][0];         vec3[1] = qed_cartesian[pos2][qed2][1];         vec3[2] = qed_cartesian[pos2][qed2][2];
  vec4[0] = grif_crystal_cartesian_110mm[ge2][0]; vec4[1] = grif_crystal_cartesian_110mm[ge2][1]; vec4[2] = grif_crystal_cartesian_110mm[ge2][2];

  cross_product(vec1,vec2,first_scattering_plane);   // the Normal vector of the plane (qed1,ge1)
  cross_product(vec3,vec4,second_scattering_plane);  // the Normal vector of the plane (qed2,ge2)
  // Now find the angle between the two scattering planes, the azimuthal
  dot = dot_product(first_scattering_plane,second_scattering_plane);
  mag = vector_magnitude_product(first_scattering_plane,second_scattering_plane);
  if(dot==0 && mag==0){ // Protection from NaN in the division
    dot = mag = 1;
  }
  angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

  // Determine the handedness based on if the scatter of the first photon is upstream or downstream
  // Downstream is positive z and positive handedness (azimuthal is positive 0 ... 180)
  // Upstream is negative z and negative handedness (azimuthal is negative -1 ... -180)
  // If the z coordinate of the HPGe is larger than z coordinate of the DSSD pixel then the scatter is in the downstream direction
  if(vec1[0]>0){
    if(vec2[2]<vec1[2]){
      //fprintf(stdout,"azimuthal_DCS, vec1 is left. Scatter is upstream\n");
      angle *= -1; }
      //  fprintf(stdout,"azimuthal_DCS, vec1 is left. Scatter is downstream\n");
    }else{
      if(vec4[2]<vec3[2]){
        //fprintf(stdout,"azimuthal_DCS, vec3 is left. Scatter is upstream\n");
        angle *= -1; }
        //fprintf(stdout,"azimuthal_DCS, vec3 is left. Scatter is downstream\n");
      }
      //printf("azimuthal_DCS angle = %f\n",angle);
      angle += 180; // angle now runs from 0-360
      return angle;
    }

    // Double Compton Scatter (DCS) - delta phi, angle between the two scattering planes
    // Calculate the azimuthal angle between the two scattering planes defined by two QED-HPGe scatter events
    // c1 and c2 are the QED pixel and HPGe of one event. These define the scattering plane.
    // c3 and c4 are the QED pixel and HPGe of one event. These define the scattering plane.
    double energy_corrected_azimuthal_DCS(int pos1, int qed1, int ge1, float ecal1, int pos2, int qed2, int ge2, float ecal2){
      double vec1[3], vec2[3], vec3[3], vec4[3], first_geometric_scattering_plane[3], second_geometric_scattering_plane[3], dot, mag, angle, scalar1, scalar2, scalar1m, scalar1p;
      double first_geometric_unit_scattering_plane[3], second_geometric_unit_scattering_plane[3], rotated_vec2[3], rotated_vec4[3], plane1[4], plane2[4], intercept1[3], intercept2[3];
      double temp_vec[3], direction1[3], direction2[3], unit1[3], unit3[3], cone_base_normal1[3],cone_base_radius1[3],tangent1[3],phi_arc1[3],o[3],cone_height1;
      double geometric_theta1, geometric_theta2, energy_derived_theta1, energy_derived_theta2, rotation_angle1, rotation_angle2;
      double enumerator1, enumerator2, denominator, boundary1[3], boundary2[3], midpoint1[3], midpoint2[3];
      double first_energy_corrected_scattering_plane[3], second_energy_corrected_scattering_plane[3];
      double cone_base_normal2[3],cone_base_radius2[3],tangent2[3],phi_arc2[3],cone_height2;

      // Calculate the scattering theta angles from the geometric centres.
      geometric_theta1 = scattering_angle_QEDGe(pos1, qed1, ge1);
      geometric_theta2 = scattering_angle_QEDGe(pos2, qed2, ge2);

      // Calculate the scattering theta angles from the Ge energy (secondary photon energy).
      // This is more accurate than the theta from geometric centres.
      energy_derived_theta1 = compton_angle(ecal1, 511.0);
      energy_derived_theta2 = compton_angle(ecal2, 511.0);

      // Adjust pos for directly accessing arrays
      pos1--; // pos1 is now 0-5 within this function
      pos2--; // pos2 is now 0-5 within this function

      // Rotate the vector of the geometric center to the energy-derived theta, keeping the same phi.
      // These are still on the same plane so would not change the delta-phi calculation.
      // Calculate where this vector intercepts the clover front-face plane to get the new coordinates which will have a different length.
      // Use this coordinate and the unit vector of the theta cone to calculate the line of interception of the clover front-face plane and theta cone conic section
      // Find the centre of this phi line with respect to the crystal boundaries.
      // Calculate the new scattering plane from this mid-point
      // Done.

      // Define the cone where the central axis is the original direction of propogation and the sides have angle energy_derived_theta1.

      // Get the plane of the Ge crystal from unit vector and central coordinate
      // grif_clover_unit_vector_110mm[(int)ge1/4];
      // Find the conic section where the Ge crystal intercepts the cone
      //

      // Calculate the vectors for the geometric scattering planes
      vec1[0] = qed_cartesian[pos1][qed1][0];         vec1[1] = qed_cartesian[pos1][qed1][1];         vec1[2] = qed_cartesian[pos1][qed1][2];
      vec2[0] = grif_crystal_cartesian_110mm[ge1][0]; vec2[1] = grif_crystal_cartesian_110mm[ge1][1]; vec2[2] = grif_crystal_cartesian_110mm[ge1][2];
      direction1[0] = vec2[0]-vec1[0]; direction1[1] = vec2[1]-vec1[1]; direction1[2] = vec2[2]-vec1[2];

      vec3[0] = qed_cartesian[pos2][qed2][0];         vec3[1] = qed_cartesian[pos2][qed2][1];         vec3[2] = qed_cartesian[pos2][qed2][2];
      vec4[0] = grif_crystal_cartesian_110mm[ge2][0]; vec4[1] = grif_crystal_cartesian_110mm[ge2][1]; vec4[2] = grif_crystal_cartesian_110mm[ge2][2];
      direction2[0] = vec3[0]-vec4[0]; direction2[1] = vec3[1]-vec4[1]; direction2[2] = vec3[2]-vec4[2];

      cross_product(vec1,vec2,first_geometric_scattering_plane);   // the Normal vector of the geometric scattering plane (qed1,ge1)
      cross_product(vec3,vec4,second_geometric_scattering_plane);  // the Normal vector of the geometric scattering plane (qed2,ge2)


      // Calculate the Normal unit vectors of the planes
      mag = sqrt( first_geometric_scattering_plane[0]*first_geometric_scattering_plane[0] + first_geometric_scattering_plane[1]*first_geometric_scattering_plane[1] + first_geometric_scattering_plane[2]*first_geometric_scattering_plane[2] );
      first_geometric_unit_scattering_plane[0] = first_geometric_scattering_plane[0] / mag;
      first_geometric_unit_scattering_plane[1] = first_geometric_scattering_plane[1] / mag;
      first_geometric_unit_scattering_plane[2] = first_geometric_scattering_plane[2] / mag;

      mag = sqrt( second_geometric_scattering_plane[0]*second_geometric_scattering_plane[0] + second_geometric_scattering_plane[1]*second_geometric_scattering_plane[1] + second_geometric_scattering_plane[2]*second_geometric_scattering_plane[2] );
      second_geometric_unit_scattering_plane[0] = second_geometric_scattering_plane[0] / mag;
      second_geometric_unit_scattering_plane[1] = second_geometric_scattering_plane[1] / mag;
      second_geometric_unit_scattering_plane[2] = second_geometric_scattering_plane[2] / mag;

      // The scattering plane does not change yet, but we need a new vector which is vec2 rotated by the rotation angle
      // The rotation axis is the Normal vector of the scattering plane
      rotation_angle1 = geometric_theta1 - energy_derived_theta1;
      rotation_angle2 = geometric_theta2 - energy_derived_theta2;

      // Rotate
      rotate_vector(direction1,first_geometric_unit_scattering_plane,rotation_angle1,rotated_vec2);
      rotate_vector(direction2,second_geometric_unit_scattering_plane,rotation_angle2,rotated_vec4);

      // Calculate the equation of the plane of the front face of this clover/crystal
      plane1[0] = grif_clover_unit_vector_110mm[(int)(ge1/4)][0];
      plane1[1] = grif_clover_unit_vector_110mm[(int)(ge1/4)][1];
      plane1[2] = grif_clover_unit_vector_110mm[(int)(ge1/4)][2];
      plane1[3] = plane1[0]*(-1*grif_crystal_cartesian_110mm[ge1][0]) + plane1[1]*(-1*grif_crystal_cartesian_110mm[ge1][1]) + plane1[2]*(-1*grif_crystal_cartesian_110mm[ge1][2]);

      // Calculate the equation of the plane of the front face of this clover/crystal
      plane2[0] = grif_clover_unit_vector_110mm[(int)(ge2/4)][0];
      plane2[1] = grif_clover_unit_vector_110mm[(int)(ge2/4)][1];
      plane2[2] = grif_clover_unit_vector_110mm[(int)(ge2/4)][2];
      plane2[3] = plane2[0]*(-1*grif_crystal_cartesian_110mm[ge2][0]) + plane2[1]*(-1*grif_crystal_cartesian_110mm[ge2][1]) + plane2[2]*(-1*grif_crystal_cartesian_110mm[ge2][2]);

      // Find the intercept point of the rotated vector and the plane of the Ge crystal face
      // Ax + By + Cz + D = 0 // plane equation of front face
      // x=o_1+d_1t, y=o_2+d_2t, z=o_3+d_3t // parametric equations of line (rotated vector)
      // A(o_1+d_1t) + B(o_2+d_2t) + C(o_3+d_3t) + D = 0 // Substitute line into plane
      // t = -1* [ A*o_1 + B*o_2 + C*o_3 + D ] / [ A*d_1 + B*d_2 + C*d_3 ] // and rearrange for t and solve.
      // Calculate intercept point using t, origin and direction vector
      scalar1 = -1* ((qed_cartesian[pos1][qed1][0]*plane1[0] + qed_cartesian[pos1][qed1][1]*plane1[1] + qed_cartesian[pos1][qed1][2]*plane1[2] + plane1[3])
      /(rotated_vec2[0]*plane1[0] + rotated_vec2[1]*plane1[1] + rotated_vec2[2]*plane1[2]) );
      intercept1[0] = qed_cartesian[pos1][qed1][0] + rotated_vec2[0]*scalar1;
      intercept1[1] = qed_cartesian[pos1][qed1][1] + rotated_vec2[1]*scalar1;
      intercept1[2] = qed_cartesian[pos1][qed1][2] + rotated_vec2[2]*scalar1;

      scalar2 = -1* ((qed_cartesian[pos2][qed2][0]*plane2[0] + qed_cartesian[pos2][qed2][1]*plane2[1] + qed_cartesian[pos2][qed2][2]*plane2[2] + plane2[3])
      /(rotated_vec4[0]*plane2[0] + rotated_vec4[1]*plane2[1] + rotated_vec4[2]*plane2[2]) );
      intercept2[0] = qed_cartesian[pos2][qed2][0] + rotated_vec4[0]*scalar2;
      intercept2[1] = qed_cartesian[pos2][qed2][1] + rotated_vec4[1]*scalar2;
      intercept2[2] = qed_cartesian[pos2][qed2][2] + rotated_vec4[2]*scalar2;

      //  Vec1 is the central axis of the cone, and Normal vector of the conic plane.
      //  Intercept is a point on the plane.
      //  Cos(theta) * (intercept - vec1) = centre of cone base plane.
      //  Vector from cone centre to intercept.
      //  Cross Product -> tangent line.

      // Define the vector of the central axis of the theta cone
      // The apex is the QED pixel coordinate
      // The height of the cone is related to the distance from the QED pixel to intercept point
      mag = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
      unit1[0] = vec1[0]/mag; unit1[1] = vec1[1]/mag; unit1[2] = vec1[2]/mag;
      cone_height1 = cos(DEGREES_TO_RADIANS*energy_derived_theta1) * sqrt(rotated_vec2[0]*rotated_vec2[0] + rotated_vec2[1]*rotated_vec2[1] + rotated_vec2[2]*rotated_vec2[2]);
      cone_base_normal1[0] = cone_height1*unit1[0]; cone_base_normal1[1] = cone_height1*unit1[1]; cone_base_normal1[2] = cone_height1*unit1[2];

      cone_base_radius1[0] = rotated_vec2[0]-cone_base_normal1[0]; cone_base_radius1[1] = rotated_vec2[1]-cone_base_normal1[1]; cone_base_radius1[2] = rotated_vec2[2]-cone_base_normal1[2];

      // The orthogonal vector will be the tangent line of the cone base at the intercept point
      cross_product(cone_base_normal1,cone_base_radius1,tangent1);
      // Convert this tangent vector to a plane and find the intercepting vector alone the Ge front face plane
      // cone_base_radius1 is the Normal vector to the tangent plane
      // intercept is a coordinate on the tangent plane
      cross_product(cone_base_radius1,(double *)grif_clover_unit_vector_110mm[(int)(ge2/4)],phi_arc1);

      // Define the vector of the central axis of the theta cone
      // The apex is the QED pixel coordinate
      // The height of the cone is related to the distance from the QED pixel to intercept point
      mag = sqrt(vec3[0]*vec3[0] + vec3[1]*vec3[1] + vec3[2]*vec3[2]);
      unit3[0] = vec3[0]/mag; unit3[1] = vec3[1]/mag; unit3[2] = vec3[2]/mag;
      cone_height2 = cos(DEGREES_TO_RADIANS*energy_derived_theta2) * sqrt(rotated_vec4[0]*rotated_vec4[0] + rotated_vec4[1]*rotated_vec4[1] + rotated_vec4[2]*rotated_vec4[2]);
      cone_base_normal2[0] = cone_height2*unit3[0]; cone_base_normal2[1] = cone_height2*unit3[1]; cone_base_normal2[2] = cone_height2*unit3[2];

      cone_base_radius2[0] = rotated_vec4[0]-cone_base_normal2[0]; cone_base_radius2[1] = rotated_vec4[1]-cone_base_normal2[1]; cone_base_radius2[2] = rotated_vec4[2]-cone_base_normal2[2];

      // The orthogonal vector will be the tangent line of the cone base at the intercept point
      cross_product(cone_base_normal2,cone_base_radius2,tangent2);
      // Convert this tangent vector to a plane and find the intercepting vector alone the Ge front face plane
      // cone_base_radius1 is the Normal vector to the tangent plane
      // intercept is a coordinate on the tangent plane
      cross_product(cone_base_radius2,(double *)grif_clover_unit_vector_110mm[(int)(ge2/4)],phi_arc2);
      //    fprintf(stdout,"phi arc2 = %.1f %.1f %.1f\n",phi_arc2[0],phi_arc2[1],phi_arc2[2]);

      // phi_arc1 is in the plane of the ge face.
      // Ge crystal is 30mm radius around the crystal centre
      // (x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - r^2 = 0 // the sphere equation
      // x=o_1+d_1t, y=o_2+d_2t, z=o_3+d_3t // parametric equations of line
      // (o_1+d_1t - x_0)^2 + (o_2+d_2t - y_0)^2 + (o_3+d_3t - z_0)^2 - r^2 = 0 // Substitute line into sphere (note the line lies within the plane)
      // (o_1+d_1*t - x_0)^2 + (o_2+d_2*t - y_0)^2 + (o_3+d_3*t - z_0)^2 - r^2 = 0
      // here o_1 is origin of line/vector, x_0 is centre of the circle
      // Treat as a Quadratic and rearrnage to find the solutions for t.
      // t = [-1*( (o_1+d_1 - x_0) + (o_2+d_2 - y_0) + (o_3+d_3 - z_0) ) +/- SQRT( ( (o_1+d_1 - x_0) + (o_2+d_2 - y_0) + (o_3+d_3 - z_0) )^2 - (d_1^2 + d_2^2 + d_3^2)*((o_1-x_0)^2+(o_2-y_0)^2+(o_3-z_0)^2-r^2) ) ]
      //     / [ (d_1^2 + d_2^2 + d_3^2 ) ]
      // Calculate intercept point using t, origin and direction vector
      //
      // First need to find an origin coordinate for the vector. The vector is phi_arc
      o[0] = rotated_vec2[0] - phi_arc1[0]*0.5; o[1] = rotated_vec2[1] - phi_arc1[1]*0.5; o[2] = rotated_vec2[2] - phi_arc1[2]*0.5;
      denominator = (phi_arc1[0]*phi_arc1[0] + phi_arc1[1]*phi_arc1[1] + phi_arc1[2]*phi_arc1[2] );
      enumerator1 = -1*( (o[0]+phi_arc1[0] - grif_crystal_cartesian_110mm[ge1][0]) + (o[1]+phi_arc1[2] - grif_crystal_cartesian_110mm[ge1][1]) + (o[2]+phi_arc1[2] - grif_crystal_cartesian_110mm[ge1][2]) );
      enumerator2 = pow( (o[0]+phi_arc1[0] - grif_crystal_cartesian_110mm[ge1][0]) + (o[1]+phi_arc1[1] - grif_crystal_cartesian_110mm[ge1][1]) + (o[2]+phi_arc1[2] - grif_crystal_cartesian_110mm[ge1][2]) ,2)
      - (phi_arc1[0]*phi_arc1[0] + phi_arc1[1]*phi_arc1[1] + phi_arc1[2]*phi_arc1[2])
      * (pow(o[0]-grif_crystal_cartesian_110mm[ge1][0],2)+pow(o[1]-grif_crystal_cartesian_110mm[ge1][1],2)+pow(o[2]-grif_crystal_cartesian_110mm[ge1][2],2) - 30*30);
      scalar1m = (enumerator1 - enumerator2) / denominator;
      scalar1p = (enumerator1 + enumerator2) / denominator;
      boundary1[0] = intercept1[0] + phi_arc1[0]*scalar1m;
      boundary1[1] = intercept1[1] + phi_arc1[1]*scalar1m;
      boundary1[2] = intercept1[2] + phi_arc1[2]*scalar1m;
      boundary2[0] = intercept1[0] + phi_arc1[0]*scalar1p;
      boundary2[1] = intercept1[1] + phi_arc1[1]*scalar1p;
      boundary2[2] = intercept1[2] + phi_arc1[2]*scalar1p;

      //  -find mid point of this line.
      midpoint1[0] = boundary1[0] + (boundary2[0]-boundary1[0])*0.5;
      midpoint1[1] = boundary1[1] + (boundary2[1]-boundary1[1])*0.5;
      midpoint1[2] = boundary1[2] + (boundary2[2]-boundary1[2])*0.5;

      o[0] = rotated_vec4[0] - phi_arc2[0]*0.5; o[1] = rotated_vec4[1] - phi_arc2[1]*0.5; o[2] = rotated_vec4[2] - phi_arc2[2]*0.5;
      denominator = (phi_arc2[0]*phi_arc2[0] + phi_arc2[1]*phi_arc2[1] + phi_arc2[2]*phi_arc2[2] );
      enumerator1 = -1*( (o[0]+phi_arc2[0] - grif_crystal_cartesian_110mm[ge2][0]) + (o[1]+phi_arc2[2] - grif_crystal_cartesian_110mm[ge2][1]) + (o[2]+phi_arc2[2] - grif_crystal_cartesian_110mm[ge2][2]) );
      enumerator2 = pow( (o[0]+phi_arc2[0] - grif_crystal_cartesian_110mm[ge2][0]) + (o[1]+phi_arc2[1] - grif_crystal_cartesian_110mm[ge2][1]) + (o[2]+phi_arc2[2] - grif_crystal_cartesian_110mm[ge2][2]) ,2)
      - (phi_arc2[0]*phi_arc2[0] + phi_arc2[1]*phi_arc2[1] + phi_arc2[2]*phi_arc2[2])
      * (pow(o[0]-grif_crystal_cartesian_110mm[ge2][0],2)+pow(o[1]-grif_crystal_cartesian_110mm[ge2][1],2)+pow(o[2]-grif_crystal_cartesian_110mm[ge2][2],2) - 30*30);
      scalar1m = (enumerator1 - enumerator2) / denominator;
      scalar1p = (enumerator1 + enumerator2) / denominator;
      boundary1[0] = intercept2[0] + phi_arc2[0]*scalar1m;
      boundary1[1] = intercept2[1] + phi_arc2[1]*scalar1m;
      boundary1[2] = intercept2[2] + phi_arc2[2]*scalar1m;
      boundary2[0] = intercept2[0] + phi_arc2[0]*scalar1p;
      boundary2[1] = intercept2[1] + phi_arc2[1]*scalar1p;
      boundary2[2] = intercept2[2] + phi_arc2[2]*scalar1p;

      //  -find mid point of this line.
      midpoint2[0] = boundary1[0] + (boundary2[0]-boundary1[0])*0.5;
      midpoint2[1] = boundary1[1] + (boundary2[1]-boundary1[1])*0.5;
      midpoint2[2] = boundary1[2] + (boundary2[2]-boundary1[2])*0.5;

      //  -Find where phi_arc1 line intercepts the boundaries of this crystal (might need to save these coordinates in a lookup table).
      //  -find mid point of this line.
      //  -This is the new intercept point.
      //  -Use this point with QED_pixel to define new vector.
      //  -Calculate new scattering plane.
      //  -Boom. Mic drop.

      // Find the centre of the phi arc of the energy-derived theta - or approximate with a straight line
      // Determine the new scattering plane from this point and the DSSD pixel
      cross_product(vec1,midpoint1,first_energy_corrected_scattering_plane);   // the Normal vector of the energy-corrected scattering plane (qed1,ge1)
      cross_product(vec3,midpoint2,second_energy_corrected_scattering_plane);  // the Normal vector of the geometric scattering plane (qed2,ge2)

      // Use this plane to calculate delta Phi
      // Now find the angle between the two scattering planes, the azimuthal
      dot = dot_product(first_geometric_scattering_plane,second_geometric_scattering_plane);
      mag = vector_magnitude_product(first_geometric_scattering_plane,second_geometric_scattering_plane);
      if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
      angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

      // Now find the angle between the two scattering planes, the azimuthal
      dot = dot_product(first_energy_corrected_scattering_plane,second_energy_corrected_scattering_plane);
      mag = vector_magnitude_product(first_energy_corrected_scattering_plane,second_energy_corrected_scattering_plane);
      if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
      angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

      // Determine the handedness based on if the scatter of the first photon is upstream or downstream
      // Downstream is positive z and positive handedness (azimuthal is positive 0 ... 180)
      // Upstream is negative z and negative handedness (azimuthal is negative -1 ... -180)
      // If the z coordinate of the HPGe is larger than z coordinate of the DSSD pixel then the scatter is in the downstream direction
      if(vec1[0]>0){
        if(vec2[2]<vec1[2]){
          //fprintf(stdout,"azimuthal_DCS, vec1 is left. Scatter is upstream\n");
          angle *= -1; }
          //  fprintf(stdout,"azimuthal_DCS, vec1 is left. Scatter is downstream\n");
        }else{
          if(vec4[2]<vec3[2]){
            //fprintf(stdout,"azimuthal_DCS, vec3 is left. Scatter is upstream\n");
            angle *= -1; }
            //fprintf(stdout,"azimuthal_DCS, vec3 is left. Scatter is downstream\n");
          }
          //printf("azimuthal_DCS angle = %f\n",angle);
          angle += 180; // angle now runs from 0-360

          if(isnan(energy_derived_theta1)||isnan(energy_derived_theta2)){
            printf("if(isNaN(energy_derived_theta1)||isNaN(energy_derived_theta2) %f,%f ==> ",energy_derived_theta1,energy_derived_theta2);
            printf("azimuthal_DCS angle = %f\n",angle);
          }

          return angle;
        }


        // Triple Compton Scatter (TCS) - delta phi, angle between the two scattering planes
        // TCS for a coincidence of [Si-Ge] and [Si-Ge-Ge] event
        // Calculate the azimuthal angle between the two scattering planes defined by a QED-HPGe plane and Ge-Ge plane
        // pos1 and ge1 are the QED pixel and HPGe of one event. These define the scattering plane.
        // The second photon undergoes an Intermediate Compton Scatter in DSSD which is ignored.
        // ge3 and ge4 are the two HPGe of the secondary Compton scatter. These define the scattering plane.
        double azimuthal_TCS_SiGe_SiGeGe(int pos1, int qed1, int ge1, int ge2, int ge3){
          double vec1[3], vec2[3], vec3[3], vec4[3], first_scattering_plane[3], second_scattering_plane[3], dot, mag, angle;

          pos1--; // pos1 is now 0-5 within this function
          vec1[0] = qed_cartesian[pos1][qed1][0];         vec1[1] = qed_cartesian[pos1][qed1][1];         vec1[2] = qed_cartesian[pos1][qed1][2];
          vec2[0] = grif_crystal_cartesian_110mm[ge1][0]; vec2[1] = grif_crystal_cartesian_110mm[ge1][1]; vec2[2] = grif_crystal_cartesian_110mm[ge1][2];

          vec3[0] = grif_crystal_cartesian_110mm[ge2][0]; vec3[1] = grif_crystal_cartesian_110mm[ge2][1]; vec3[2] = grif_crystal_cartesian_110mm[ge2][2];
          vec4[0] = grif_crystal_cartesian_110mm[ge3][0]; vec4[1] = grif_crystal_cartesian_110mm[ge3][1]; vec4[2] = grif_crystal_cartesian_110mm[ge3][2];

          cross_product(vec1,vec2,first_scattering_plane);   // the Normal vector of the plane (qed1,ge1)
          cross_product(vec3,vec4,second_scattering_plane);  // the Normal vector of the plane (ge2,ge3)
          // Now find the angle between the two scattering planes, the azimuthal
          dot = dot_product(first_scattering_plane,second_scattering_plane);
          mag = vector_magnitude_product(first_scattering_plane,second_scattering_plane);
          if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
          angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

          // Determine the handedness based on if the scatter of the first photon is upstream or downstream
          // Downstream is positive z and positive handedness (azimuthal is positive 0 ... 180)
          // Upstream is negative z and negative handedness (azimuthal is negative -1 ... -180)
          // If the z coordinate of the HPGe is larger than z coordinate of the DSSD pixel then the scatter is in the downstream direction
          if(vec1[0]>0){
            if(vec2[2]<vec1[2]){ angle *= -1; }
          }else{
            if(vec4[2]<vec3[2]){ angle *= -1; }
          }
          //printf("azimuthal_DCS angle = %f\n",angle);
          angle += 180; // angle now runs from 0-360
          return angle;
        }

        // Double Compton Scatter (TCS) - delta phi, angle between the two scattering planes
        // Calculate the azimuthal angle between the two scattering planes defined by two QED-HPGe scatter events
        // c1 and c2 are the QED pixel and HPGe of one event. These define the scattering plane.
        // The second photon undergoes an Intermediate Compton Scatter in DSSD which is ignored.
        // c3 and c4 are the two HPGe of the secondary Compton scatter. These define the scattering plane.
        double azimuthal_TCS_GeGe_SiGeGe(int ge1, int ge2, int ge3, int ge4){
          double vec1[3], vec2[3], vec3[3], vec4[3], first_scattering_plane[3], second_scattering_plane[3], dot, mag, angle;

          vec1[0] = grif_crystal_cartesian_110mm[ge1][0]; vec1[1] = grif_crystal_cartesian_110mm[ge1][1]; vec1[2] = grif_crystal_cartesian_110mm[ge1][2];
          vec2[0] = grif_crystal_cartesian_110mm[ge2][0]; vec2[1] = grif_crystal_cartesian_110mm[ge2][1]; vec2[2] = grif_crystal_cartesian_110mm[ge2][2];

          vec3[0] = grif_crystal_cartesian_110mm[ge3][0]; vec3[1] = grif_crystal_cartesian_110mm[ge3][1]; vec3[2] = grif_crystal_cartesian_110mm[ge3][2];
          vec4[0] = grif_crystal_cartesian_110mm[ge4][0]; vec4[1] = grif_crystal_cartesian_110mm[ge4][1]; vec4[2] = grif_crystal_cartesian_110mm[ge4][2];

          cross_product(vec1,vec2,first_scattering_plane);   // the Normal vector of the plane (qed1,ge1)
          cross_product(vec3,vec4,second_scattering_plane);  // the Normal vector of the plane (qed2,ge2)
          // Now find the angle between the two scattering planes, the azimuthal
          dot = dot_product(first_scattering_plane,second_scattering_plane);
          mag = vector_magnitude_product(first_scattering_plane,second_scattering_plane);
          if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
          angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

          // Determine the handedness based on if the scatter of the first photon is upstream or downstream
          // Downstream is positive z and positive handedness (azimuthal is positive 0 ... 180)
          // Upstream is negative z and negative handedness (azimuthal is negative -1 ... -180)
          // If the z coordinate of the HPGe is larger than z coordinate of the DSSD pixel then the scatter is in the downstream direction
          if(vec1[0]>0){
            if(vec2[2]<vec1[2]){ angle *= -1; }
          }else{
            if(vec4[2]<vec3[2]){ angle *= -1; }
          }
          angle += 180; // angle now runs from 0-360
          return angle;
        }

        // Triple Compton Scatter (TCS) - delta phi, angle between the two scattering planes
        // TCS for a coincidence of [Ge-Ge] and [Si-Si-Ge] event
        // Calculate the azimuthal angle between the two scattering planes defined by a QED-HPGe plane and HPGe-HPGe plane
        // ge2 and ge3 are the HPGe of one Compton event. These define the scattering plane.
        // The second photon undergoes an Intermediate Compton Scatter in DSSD which is ignored.
        // pos1,qed1 and ge1 are the QED pixel and HPGe of the secondary Compton scatter. These define the scattering plane.
        double azimuthal_TCS_GeGe_SiSiGe(int pos1, int qed1, int ge1, int ge2, int ge3){
          double vec1[3], vec2[3], vec3[3], vec4[3], first_scattering_plane[3], second_scattering_plane[3], dot, mag, angle;

          pos1--; // pos are now 0-5 within this function
          vec1[0] = qed_cartesian[pos1][qed1][0];         vec1[1] = qed_cartesian[pos1][qed1][1];         vec1[2] = qed_cartesian[pos1][qed1][2];
          vec2[0] = grif_crystal_cartesian_110mm[ge1][0]; vec2[1] = grif_crystal_cartesian_110mm[ge1][1]; vec2[2] = grif_crystal_cartesian_110mm[ge1][2];

          vec3[0] = grif_crystal_cartesian_110mm[ge2][0]; vec3[1] = grif_crystal_cartesian_110mm[ge2][1]; vec3[2] = grif_crystal_cartesian_110mm[ge2][2];
          vec4[0] = grif_crystal_cartesian_110mm[ge3][0]; vec4[1] = grif_crystal_cartesian_110mm[ge3][1]; vec4[2] = grif_crystal_cartesian_110mm[ge3][2];

          cross_product(vec1,vec2,first_scattering_plane);   // the Normal vector of the plane (qed1,ge1)
          cross_product(vec3,vec4,second_scattering_plane);  // the Normal vector of the plane (qed3,ge2)
          // Now find the angle between the two scattering planes, the azimuthal
          dot = dot_product(first_scattering_plane,second_scattering_plane);
          mag = vector_magnitude_product(first_scattering_plane,second_scattering_plane);
          if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
          angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

          // Determine the handedness based on if the scatter of the first photon is upstream or downstream
          // Downstream is positive z and positive handedness (azimuthal is positive 0 ... 180)
          // Upstream is negative z and negative handedness (azimuthal is negative -1 ... -180)
          // If the z coordinate of the HPGe is larger than z coordinate of the DSSD pixel then the scatter is in the downstream direction
          if(vec1[0]>0){
            if(vec2[2]<vec1[2]){ angle *= -1; }
          }else{
            if(vec4[2]<vec3[2]){ angle *= -1; }
          }
          angle += 180; // angle now runs from 0-360
          return angle;
        }

        // Triple Compton Scatter (TCS) - delta phi, angle between the two scattering planes
        // TCS for a coincidence of [Si-Ge] and [Si-Si-Ge] event
        // Calculate the azimuthal angle between the two scattering planes defined by a QED-HPGe plane and QED-HPGe plane
        // pos1,qed1 and ge1 are the QED pixel and HPGe of one Compton event. These define the scattering plane.
        // The second photon undergoes an Intermediate Compton Scatter in DSSD (pos2,qed2) which is ignored.
        // pos3,qed3 and ge2 are the QED pixel and HPGe of the secondary Compton scatter. These define the scattering plane.
        double azimuthal_TCS_SiGe_SiSiGe(int pos1, int qed1, int ge1, int pos3, int qed3, int ge2){
          double vec1[3], vec2[3], vec3[3], vec4[3], first_scattering_plane[3], second_scattering_plane[3], dot, mag, angle;

          pos1--; pos3--; // pos are now 0-5 within this function
          vec1[0] = qed_cartesian[pos1][qed1][0];         vec1[1] = qed_cartesian[pos1][qed1][1];         vec1[2] = qed_cartesian[pos1][qed1][2];
          vec2[0] = grif_crystal_cartesian_110mm[ge1][0]; vec2[1] = grif_crystal_cartesian_110mm[ge1][1]; vec2[2] = grif_crystal_cartesian_110mm[ge1][2];

          vec3[0] = qed_cartesian[pos3][qed3][0];         vec3[1] = qed_cartesian[pos3][qed3][1];         vec3[2] = qed_cartesian[pos3][qed3][2];
          vec4[0] = grif_crystal_cartesian_110mm[ge2][0]; vec4[1] = grif_crystal_cartesian_110mm[ge2][1]; vec4[2] = grif_crystal_cartesian_110mm[ge2][2];

          cross_product(vec1,vec2,first_scattering_plane);   // the Normal vector of the plane (qed1,ge1)
          cross_product(vec3,vec4,second_scattering_plane);  // the Normal vector of the plane (qed3,ge2)
          // Now find the angle between the two scattering planes, the azimuthal
          dot = dot_product(first_scattering_plane,second_scattering_plane);
          mag = vector_magnitude_product(first_scattering_plane,second_scattering_plane);
          if(dot==0 && mag==0){ dot = mag = 1; } // Protection from NaN in the division
          angle = RADIANS_TO_DEGREES*fast_acos( dot / mag ); // acos is very slow

          // Determine the handedness based on if the scatter of the first photon is upstream or downstream
          // Downstream is positive z and positive handedness (azimuthal is positive 0 ... 180)
          // Upstream is negative z and negative handedness (azimuthal is negative -1 ... -180)
          // If the z coordinate of the HPGe is larger than z coordinate of the DSSD pixel then the scatter is in the downstream direction
          if(vec1[0]>0){
            if(vec2[2]<vec1[2]){ angle *= -1; }
          }else{
            if(vec4[2]<vec3[2]){ angle *= -1; }
          }
          angle += 180; // angle now runs from 0-360
          return angle;
        }
