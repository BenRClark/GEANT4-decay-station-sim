#ifndef DetectionSystemTigress_h
#define DetectionSystemTigress_h 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "globals.hh"
#include "G4SDManager.hh"
//#include "DetectableSD.hh"
//#include "EventDataOutput.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4AssemblyVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "Materials.hh"

///////////////////////////////////////////////////////////////////////
// DetectionSystemTigress
///////////////////////////////////////////////////////////////////////
class DetectionSystemTigress
{
  public:
  DetectionSystemTigress(G4SDManager* SDman, /*EventDataOutput* edo,*/ G4int numberOfDetectors, G4int detectorPosition, G4bool hevimet);
    ~DetectionSystemTigress();

  public:
    void Build(G4LogicalVolume* exp_hall_log);  
  //    DetectableSD* getGermaniumSD(){return germanium_block_SD;};      
  //DetectableSD* getLeftCasingSD(){return left_casing_SD;};  
  // DetectableSD* getRightCasingSD(){return right_casing_SD;};  
  // DetectableSD* getLeftExSD(){return left_extension_SD;};  
  // DetectableSD* getRightExSD(){return right_extension_SD;};  
  // DetectableSD* getBackPlugSD(){return back_plug_SD;};         

  private:    
  //  DetectableSD* crystal_block_SD;
    G4LogicalVolume* expHallLog; 
    G4SDManager* mySDman;
  //    EventDataOutput* myEdo;
    
    G4RotationMatrix* rotate_null;
    G4ThreeVector move_null;

    double coords[20][5];
    
    G4String suppressor_shape;
    G4String hevimet_choice;

    //Settings for the suppressor shield
    G4double thickness_double;
    G4double side_thickness_double;
    G4double extension_thickness_double;
    G4String back_material;

    G4double radial_distance; 
    G4int number_of_detectors;
    G4int number_of_segments;
    G4double origin_to_crystal_distance;
    G4bool dead_layer_include_flag;
    G4double inner_dead_layer_thickness;
    G4double outer_dead_layer_thickness;
 
    // Suppressor design flags
    G4bool suppressor_shells_include_flag;
   
    //Cold Finger
    G4double cold_finger_outer_shell_radius;
    G4double cold_finger_shell_thickness;
    G4double cold_finger_shell_length;	
    G4double rear_plate_thickness;
    G4double cold_finger_radius;
    G4double cold_finger_space;
    G4double cold_finger_length;
    G4double coolant_length;
    G4double coolant_radius;
    
    G4double fet_air_hole_radius;
    G4double aluminum_cold_finger_thickness;
    G4double aluminum_cold_finger_radius;
    
    G4int suppressor_design;
    G4int suppressor_position_selector;
    G4int hevimet_selector;  

    G4double forward_inner_radius;
    G4double back_inner_radius;

    // For the optimization of the depth segmentation
    G4double depth_segmentation_adjustment;
    
    //the germanium detector's variables (for one "leaf" of the clover)
    G4double germanium_outer_radius;
    G4double germanium_hole_radius;
    G4double germanium_width;
    G4double germanium_length;
    G4double germanium_hole_dist_from_face;
    G4double germanium_dist_from_can_face;
    G4double germanium_bent_length;
    G4double germanium_shift;    		//the amount by which the two sides adjacent to the other 
				 		//germanium crystals are cut closer to the center
				 	
    G4double germanium_separation;  		//the space between the quarter detectors    
    G4double germanium_corner_cone_end_length; 	//the ending length of the cones
                                               	//at the corners of each leaf
    
    //basic germanium detector's values
    G4double detector_block_length;
    G4double detector_block_height;
    G4double detector_block_trapezoidal_height;

    G4double detector_total_width;
    G4double detector_total_length;
    G4double bent_end_angle;
    G4double bent_end_length;
    G4double can_face_thickness;
    G4double can_side_thickness;

    G4double hevimet_tip_thickness;
    G4double hevimet_tip_angle;

    //Values for the BGO
    G4double suppressor_cut_extension;
    G4double suppressor_shell_thickness;
    
    G4double back_BGO_thickness;
    G4double BGO_chopped_tip;
    G4double BGO_movable_space;
    G4double side_suppressor_back_offset;
    G4double side_BGO_thickness;
    G4double BGO_can_seperation;
    G4double side_BGO_length;
    G4double side_suppressor_length;
    G4double BGO_trap_length;
    G4double suppressor_extension_thickness;
    G4double suppressor_extension_length;
    G4double suppressor_extension_angle;

    G4String back_suppressor_material;

    //Values for the HeavyMet
    G4double HeavyMet_thickness;
    G4double HeavyMet_inside_angle;
    

    G4double air_box_front_width;
    G4double air_box_front_length;
    G4double air_box_back_length;

    G4double shift;

    G4int germanium_copy_number;   
    G4int left_suppressor_side_copy_number;
    G4int right_suppressor_side_copy_number;
    G4int left_suppressor_extension_copy_number;
    G4int right_suppressor_extension_copy_number;
    G4int back_suppressor_copy_number; 
    
    G4double rhombi_diameter;
    G4double new_rhombi_radius;
    G4double detector_position_shift;
    G4double applied_back_shift;

    G4int germanium_selector;
    G4int can_selector;
    G4int BGO_selector;
    G4int cold_finger_selector;
   
  private:    
  //DetectableSD* germanium_block_SD;
  // DetectableSD* left_casing_SD;
  // DetectableSD* right_casing_SD;
  // DetectableSD* left_extension_SD;
  // DetectableSD* right_extension_SD;
  // DetectableSD* back_plug_SD;

  private:
    // Assembly volumes
    G4AssemblyVolume* assembly;
    G4AssemblyVolume* germaniumAssembly;
    G4AssemblyVolume* leftSuppressorCasingAssembly;
    G4AssemblyVolume* rightSuppressorCasingAssembly;
    G4AssemblyVolume* leftSuppressorExtensionAssembly;
    G4AssemblyVolume* rightSuppressorExtensionAssembly;
    G4AssemblyVolume* suppressorBackAssembly;
      
    // Logical volumes
    G4LogicalVolume* air_box_log;
 
    // methods to construct all of the components of the detector
    void ConstructNewSuppressorCasingWithShells();
    void ConstructComplexDetectorBlockWithDeadLayer();
    void ConstructDetector();
    void ConstructComplexDetectorBlock();
    void ConstructColdFinger();
    void ConstructNewSuppressorCasing();
    void ConstructNewHeavyMet();

  private:

    //LogicalVolumes used in ConstructBasicDetectorBlock
    G4LogicalVolume* germanium_block_log;

    //Logical Volumes used in ConstructComplexDetectorBlock:
    G4LogicalVolume* germanium_block1_log;
    G4LogicalVolume* germanium_hole_log;
    G4LogicalVolume* inner_dead_layer_log;
    G4LogicalVolume* inner_dead_layer_cap_log;
    G4LogicalVolume* outer_dead_layer_log;
    
    //Logical Volumes used in ConstructBGOCasing:
    G4LogicalVolume* back_BGO_log;
    G4LogicalVolume* BGO_casing_log;

    //Logical Volumes used in ConstructNewSuppressorCasing:
    G4LogicalVolume* shell_for_back_quarter_suppressor_log;
    G4LogicalVolume* shell_for_right_suppressor_log;
    G4LogicalVolume* shell_for_left_suppressor_log;
    G4LogicalVolume* shell_for_right_suppressor_extension_log;
    G4LogicalVolume* shell_for_left_suppressor_extension_log;

    G4LogicalVolume* cap_for_right_suppressor_log;

    G4LogicalVolume* back_quarter_suppressor_log;
    G4LogicalVolume* right_suppressor_log;
    G4LogicalVolume* left_suppressor_log;
    G4LogicalVolume* right_suppressor_extension_log;
    G4LogicalVolume* left_suppressor_extension_log;
    
  private:
    //Logical Volumes used in ConstructDetector:
    G4LogicalVolume* front_face_log;
    G4LogicalVolume* right_bent_piece_log;
    G4LogicalVolume* left_bent_piece_log;
    G4LogicalVolume* top_bent_piece_log;
    G4LogicalVolume* bottom_bent_piece_log;
    G4LogicalVolume* right_wedge_log;
    G4LogicalVolume* left_wedge_log;
    G4LogicalVolume* top_wedge_log;
    G4LogicalVolume* bottom_wedge_log;
    G4LogicalVolume* upper_right_cone_log;
    G4LogicalVolume* lower_right_cone_log;
    G4LogicalVolume* upper_left_cone_log;
    G4LogicalVolume* lower_left_cone_log;
    G4LogicalVolume* upper_right_tube_log;
    G4LogicalVolume* lower_right_tube_log;
    G4LogicalVolume* upper_left_tube_log;
    G4LogicalVolume* lower_left_tube_log;
    G4LogicalVolume* right_side_panel_log;
    G4LogicalVolume* left_side_panel_log;
    G4LogicalVolume* top_side_panel_log;
    G4LogicalVolume* bottom_side_panel_log;
    G4LogicalVolume* rear_plate_log;
    G4LogicalVolume* finger_shell_log;
    G4LogicalVolume* tank_log;

    //Logical Volumes used in ConstructColdFinger:
    G4LogicalVolume* finger_log;
    G4LogicalVolume* extra_cold_block_log;
    G4LogicalVolume* fet_air_hole_log;
    G4LogicalVolume* aluminum_cold_finger_log;
    
    //Logical Volumes used in ConstructNewHeavyMet:
    G4LogicalVolume* hevimet_log;

    //internal methods for ConstructCan()
    G4Box* squareFrontFace();
    G4Trap* cornerWedge();
    G4Para* bentSidePiece();
    G4Box* otherBentSidePiece();
    G4Cons* roundedEndEdge();
    G4Tubs* cornerTube();
    G4Box* sidePanel();
    G4SubtractionSolid* rearPlate();
    G4Tubs* coldFingerShell();
    G4Tubs* liquidNitrogenTank();

    //internal methods for ConstructBasicDetectorBlock()
    G4Box* rectangularSegment();
    G4Trd* trapezoidalSegment();

    //internal methods for ConstructComplexDetectorBlock()
    G4SubtractionSolid* quarterDetector(); 

    //internal methods for ConstructColdFinger()
    G4Tubs* airHole();
    G4Box* endPlate();
    G4Tubs* finger();
    G4Tubs* aluminumColdFinger();
    
    //internal methods for ConstructBGOCasing()
    G4SubtractionSolid* backBGO();
    G4SubtractionSolid* BGOCasing();
    G4SubtractionSolid* frontBGO();
    G4Trd* sideBGO();

    //internal methods for ConstructNewSuppressorCasing()
    G4SubtractionSolid* backSuppressorQuarter();
    G4SubtractionSolid* frontRightSlantSuppressor();
    G4SubtractionSolid* frontLeftSlantSuppressor();
    G4SubtractionSolid* rightSuppressorExtension();
    G4SubtractionSolid* leftSuppressorExtension();

    G4SubtractionSolid* choppingFrontRightSlantSuppressor();
    G4SubtractionSolid* choppingFrontLeftSlantSuppressor();
    G4SubtractionSolid* choppingRightSuppressorExtension();
    G4SubtractionSolid* choppingLeftSuppressorExtension();

    //internal methods for New SuppressorCasingWithShells
    G4SubtractionSolid* shellForBackSuppressorQuarter();
    G4SubtractionSolid* shellForFrontRightSlantSuppressor();
    G4SubtractionSolid* shellForFrontLeftSlantSuppressor();
    G4SubtractionSolid* shellForRightSuppressorExtension();
    G4SubtractionSolid* shellForLeftSuppressorExtension();

    //internal methods for ConstructNewHeavyMet()
    G4SubtractionSolid* newHeavyMet();

           
  private: 
    G4String crystal_material;
    G4String can_material;
    G4String vacuum_material;
    G4double detail_view_end_angle;
    G4double crystal_length_x;   
    G4double crystal_length_y;
    G4double crystal_length_z;
    G4double crystal_inner_radius;
    G4double crystal_outer_radius;
    G4double can_thickness;
    G4double can_inner_radius;
    G4double can_lid_inner_radius;
    G4double can_lid_outer_radius;
    G4double can_front_lid_thickness;
    G4double can_back_lid_thickness;
    G4double can_face_dist_from_origin;
    G4double crystal_dist_from_can_face;
    G4double crystal_dist_from_can_back;
    G4double can_length_z;
    G4double crystal_dist_from_origin;
        
  private: 
    // internal methods
    void BuildOneDetector();  
    void PlaceDetector(G4int detector_number);
        
  private:
    G4ThreeVector GetDirection(G4double theta, G4double phi);
};

#endif

