

# Uploading large files to GitHub commands

cd "OneDrive - The University of Queensland\Desktop\GitHub"

git lfs install
# Should return Git LFS intialized

git lfs track "*.tif"
git lfs track "*.tar"
git lfs track "*.xml"

git lfs push --all origin main # Only works for the first file, skip this line for every subsequent file

git add .

git commit -m "Upload large sentinel fire data file"

git push -u origin main