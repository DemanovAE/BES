/**
 * \class StFemtoArrays
 * \brief Keeps FemtoDst names and types
 *
 * The class holds FemtoDst names and types that
 * are used for branches and other things
 *
 * \author Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date July 15, 2018
 */

#ifndef StFemtoArrays_h
#define StFemtoArrays_h

//_________________
class StFemtoArrays {

 public:
  /// Default constructor
  StFemtoArrays();
  /// Destructor
  ~StFemtoArrays();

  /// Number of used arrays
  enum {NAllFemtoArrays = 5};

  /// Names of the TBranches in the TTree/File
  static const char* femtoArrayNames[NAllFemtoArrays];

  /// Names of the classes, the TClonesArrays are arrays of this type
  static const char* femtoArrayTypes[NAllFemtoArrays];

  /// Maximum sizes of the TClonesArrays
  static int femtoArraySizes[NAllFemtoArrays];

  /// Setup type indices
  enum TypeIndex { Event=0, Track, V0, Xi, EpdHit };
};

#endif //StFemtoArrays_h
