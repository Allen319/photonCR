#include <FileInPath.h>

#include <cstdlib>
#include <iostream>
#include <boost/algorithm/string/predicate.hpp>


namespace fs = std::filesystem;


void FileInPath::AddLocation(fs::path const &path) {
  GetInstance().locations_.emplace_back(path);
}


fs::path FileInPath::Resolve(fs::path const &subDir, fs::path const &path) {

  // Check if an absolute or explicit relative path is provided. In this case
  // the first argument is ignored and the path is returned as is.
  auto const pathString = path.string();

  if (boost::starts_with(pathString, "/") or
      boost::starts_with(pathString, "./") or
      boost::starts_with(pathString, "../")) {

    // Make sure the requested file exists

    return path;
  }


  // Loop over all possible locations, giving preference to ones added later
  auto const &locations = GetInstance().locations_;

  for (auto locationIt = locations.crbegin(); locationIt != locations.crend();
      ++locationIt) {
    fs::path tryPath;

    // Try to resolve the path using the provided subdirectory
    if (not subDir.empty()) {
      tryPath = *locationIt / subDir / path;

      if (fs::exists(tryPath))
        return tryPath;
    }

    // Try to resolve the path ignoring the subdir
    tryPath = *locationIt / path;

    if (fs::exists(tryPath))
      return tryPath;
  }


  // If none of the above succeeds, try to resolve the path with respect to the
  // currect working directory, with and without the subdirectory
  if (not subDir.empty()) {
    fs::path tryPath = subDir / path;

    if (fs::exists(tryPath))
      return fs::current_path() / tryPath;
  }

  if (fs::exists(path))
    return fs::current_path() / path;

  std::cout << "error"<< std::endl;
  // If the workflow has reached this point, the path has not been resolved
}


fs::path FileInPath::Resolve(fs::path const &path) {
  return Resolve({}, path);
}

int FileInPath::GetFilenames(std::filesystem::path const &path,
                          std::vector<std::string> & filenames,
                          bool recur,
                          std::string suffix) {
  if (! fs::exists(path))
  {
    return -1;
  }
  fs::directory_iterator end_iter;
  for (fs::directory_iterator iter(path); iter != end_iter; ++iter) {
    if (fs::is_regular_file(iter->status()))
    {
      if (suffix == "")
        filenames.push_back(iter->path().string());
      else if (boost::algorithm::ends_with(iter->path().string(), "." + suffix))
        filenames.push_back(iter->path().string());
    }
    if (recur) {
      if (fs::is_directory(iter->status()))
        GetFilenames(iter->path().string(), filenames, recur, suffix);
    }
  }
  return filenames.size();
}

FileInPath::FileInPath() {
  // Read the install path from the environment
  char const *installPath = std::getenv("CMAKE_BINARY_DIR");

  // Specify default locations. Cannot use method AddLocation here since it will
  // indirectly call the constructor.
  locations_.emplace_back(fs::path{installPath} / "config");
  locations_.emplace_back(fs::path{installPath} / "samples");
  locations_.emplace_back(fs::path{installPath} / "environment");
}


FileInPath &FileInPath::GetInstance() {
  static FileInPath instance;
  return instance;
}
