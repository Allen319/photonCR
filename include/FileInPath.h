#ifndef FILEINPATH_H_
#define FILEINPATH_H_

#include <filesystem>
#include <vector>


/**
 * \brief Allows to resolve a (possibly) relative path with respect to several
 * potential locations
 *
 * A file or directory is search for in several locations, similarly to how
 * \c $PATH is used. The following locations are included by default:
 *  - \c $HZZ2L2NU_BASE/config
 *  - \c $HZZ2L2NU_BASE/data
 *
 * Other locations can be added with \ref AddLocation.
 *
 * If a path does not correspond to an existing file or directory, an exception
 * is thrown.
 * 
 * This class is a singleton, and user cannot construct an instance of it.
 * Instead, all functionality is implemented in static methods.
 */
class FileInPath {
 public:
  /**
   * \brief Adds a new location in which files with be searched
   * 
   * The new location takes preference over all paths added previously.
   */
  static void AddLocation(std::filesystem::path const &path);
  
  /**
   * \brief Resolves a path, allowing for an optional subdirectory
   * 
   * If the path starts with "/", "./", or "../", it is treated as an absolute
   * or explicit relative path and returned unchanched after verifying that such
   * file or directory exists. Otherwise the method tries to resolve it with
   * respect to all defined locations, in a reversed order of their definition.
   * For each location, provided subdirectory is first added to it, and the
   * resolution is attempted. If such file or directory is not found, the
   * subdirectory is omitted, and the resolution is attempted again. Finally,
   * the path is searched for in the current working directory (the one in which
   * the executable is being run), with and without the subdirectory. If all
   * attempts to find the path fail, an exception is thrown.
   */
  static std::filesystem::path Resolve(std::filesystem::path const &subDir,
                                       std::filesystem::path const &path);
  
  /**
   * \brief Resolves a path
   * 
   * Works in the same way as the other version but does not include the
   * additional subdirectory.
   */
  static std::filesystem::path Resolve(std::filesystem::path const &path);

  /**
   * \brief Get paths of files under a cerntain path
   * 
   * Get all absolute paths of files under a dir. One can choose include the
   * additional subdirectory or not, as well as the suffix of the files.
   */
  static int GetFilenames(std::filesystem::path const &dir, 
                          std::vector<std::string> & filenames,
                          bool recur = false, std::string suffix = "");
 private:
  /**
   * \brief Constructor
   * 
   * The constructor is private because the class is a singletop. It reads the
   * value of environmental variable \c HZZ2L2NU_BASE and sets default
   * locations. If the variable is not set, an exception is thrown.
   */
  FileInPath();

  FileInPath(FileInPath const &) = delete;
  FileInPath &operator=(FileInPath const &) = delete;

  /// Returns the only instance of this singleton
  static FileInPath &GetInstance();

  /// Locations with respect to which paths are resolved
  std::vector<std::filesystem::path> locations_;
};

#endif  // FILEINPATH_H_

