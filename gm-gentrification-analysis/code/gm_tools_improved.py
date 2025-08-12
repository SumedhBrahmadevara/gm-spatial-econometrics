#!/usr/bin/env python3
"""
gm_tools.py — CLI toolkit for Greater Manchester spatial data preparation and analysis.

This module provides utilities for processing spatial data related to Greater Manchester,
including OpenStreetMap feature extraction, spatial control variable computation,
and distance calculations.

Subcommands:
  extract-osm                  Download OSM features for a bounding box (greenspace/schools/stations)
  spatial-controls             Build distance to station, percentage greenspace, and school count metrics
  area                         Compute MSOA11 area in square kilometers
  dist-cbd                     Calculate distance from MSOA centroid to central business district
  dist-station-over-time       Calculate distance to nearest station for multiple years

Expected directory structure (customizable via CLI flags):
  data/
    raw/
      shapefiles/infuse_msoa_lyr_2011.shp
      osm/ (OSM outputs directory)
    reference/
      msoa_codes.txt                    # one MSOA11 code per line
      station_opening_years.csv         # name,opening_year
    processed/
      spatial/

Dependencies:
  - pandas
  - geopandas
  - osmnx (optional, required only for extract-osm command)
  - shapely
  - pyproj

Author: [Your Name]
Date: [Date]
Version: 1.0.0
"""

from __future__ import annotations
from pathlib import Path
import argparse
import logging
import sys
from typing import List, Set, Tuple

import pandas as pd
import geopandas as gpd

# Optional dependency for OSM extraction
try:
    import osmnx as ox
    OSM_AVAILABLE = True
except ImportError:
    ox = None
    OSM_AVAILABLE = False

# Coordinate reference system constants
EPSG_WGS84 = 4326      # WGS84 (latitude/longitude)
EPSG_BNG = 27700       # British National Grid (meters)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class SpatialDataError(Exception):
    """Custom exception for spatial data processing errors."""
    pass


def get_base_directory() -> Path:
    """
    Get the base directory of the project (parent of script directory).
    
    Returns:
        Path: Base directory path
    """
    return Path(__file__).resolve().parents[1]


def read_msoa_codes(file_path: Path) -> List[str]:
    """
    Read MSOA codes from a text file.
    
    Args:
        file_path: Path to file containing MSOA codes (one per line)
        
    Returns:
        List of MSOA codes
        
    Raises:
        FileNotFoundError: If the codes file doesn't exist
        SpatialDataError: If the file is empty or invalid
    """
    if not file_path.exists():
        raise FileNotFoundError(f"MSOA codes file not found: {file_path}")
    
    try:
        codes = [line.strip() for line in file_path.read_text().splitlines() if line.strip()]
        if not codes:
            raise SpatialDataError(f"No valid MSOA codes found in {file_path}")
        logger.info(f"Loaded {len(codes)} MSOA codes from {file_path}")
        return codes
    except Exception as e:
        raise SpatialDataError(f"Error reading MSOA codes from {file_path}: {e}")


def validate_required_columns(gdf: gpd.GeoDataFrame, required_columns: List[str]) -> None:
    """
    Validate that a GeoDataFrame contains required columns.
    
    Args:
        gdf: GeoDataFrame to validate
        required_columns: List of required column names
        
    Raises:
        SpatialDataError: If required columns are missing
    """
    missing_columns = [col for col in required_columns if col not in gdf.columns]
    if missing_columns:
        raise SpatialDataError(f"Missing required columns: {missing_columns}")


def extract_osm_features(args: argparse.Namespace) -> None:
    """
    Extract OpenStreetMap features for specified bounding box.
    
    Args:
        args: Command line arguments containing bbox and output directory
        
    Raises:
        RuntimeError: If osmnx is not available
        SpatialDataError: If OSM extraction fails
    """
    if not OSM_AVAILABLE:
        raise RuntimeError(
            "osmnx library is required for extract-osm command. "
            "Install with: pip install osmnx"
        )

    output_dir: Path = args.outdir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Bounding box: (north, south, east, west)
    bbox = (args.north, args.south, args.east, args.west)
    logger.info(f"Extracting OSM features for bbox: {bbox}")

    # Define OSM feature tags
    feature_definitions = {
        "greenspaces.geojson": {
            "leisure": ["park", "recreation_ground", "pitch"],
            "landuse": ["grass", "meadow"],
            "natural": ["wood"]
        },
        "schools.geojson": {
            "amenity": ["school", "college", "university"]
        },
        "transport_stations.geojson": {
            "railway": ["station", "tram_stop"]
        }
    }

    for filename, tags in feature_definitions.items():
        try:
            logger.info(f"Extracting {filename}...")
            gdf = ox.features_from_bbox(*bbox, tags=tags)
            output_path = output_dir / filename
            gdf.to_file(output_path, driver="GeoJSON")
            logger.info(f"Successfully saved {filename} with {len(gdf)} features")
        except Exception as e:
            logger.error(f"Failed to extract {filename}: {e}")
            raise SpatialDataError(f"OSM extraction failed for {filename}: {e}")


def compute_spatial_controls(args: argparse.Namespace) -> None:
    """
    Compute spatial control variables for MSOAs.
    
    Creates the following variables:
    - dist_to_station_km: Distance to nearest transport station
    - pct_greenspace: Percentage of area that is greenspace
    - school_count_1km: Number of schools within 1km radius
    
    Args:
        args: Command line arguments with input/output paths
    """
    args.out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Computing spatial control variables...")

    try:
        # Load and prepare MSOA data
        logger.info("Loading MSOA boundaries...")
        msoas = gpd.read_file(args.shapefile).to_crs(epsg=EPSG_BNG)
        valid_codes = set(read_msoa_codes(args.gm_codes))
        
        gm_msoas = msoas[msoas["geo_code"].isin(valid_codes)].copy()
        validate_required_columns(gm_msoas, ["geo_code", "geometry"])
        
        gm_msoas = gm_msoas.rename(columns={"geo_code": "msoa11cd"})
        gm_msoas["centroid"] = gm_msoas.geometry.centroid
        logger.info(f"Processing {len(gm_msoas)} MSOA areas")

        # Load transport stations
        logger.info("Loading transport stations...")
        stations = gpd.read_file(args.stations).to_crs(epsg=EPSG_BNG)
        stations["geometry"] = stations.geometry.centroid
        logger.info(f"Loaded {len(stations)} stations")

        # Load greenspace data
        logger.info("Loading greenspace data...")
        greenspaces = gpd.read_file(args.greenspace).to_crs(epsg=EPSG_BNG)
        logger.info(f"Loaded {len(greenspaces)} greenspace features")

        # Load school data
        logger.info("Loading school data...")
        schools = gpd.read_file(args.schools).to_crs(epsg=EPSG_BNG)
        schools["geometry"] = schools.geometry.centroid
        logger.info(f"Loaded {len(schools)} schools")

        # Calculate distance to nearest station
        logger.info("Calculating distances to nearest stations...")
        nearest_stations = gpd.sjoin_nearest(
            gm_msoas[["msoa11cd", "geometry"]],
            stations[["geometry"]],
            how="left",
            distance_col="dist_to_station_m",
        )
        gm_msoas = gm_msoas.merge(
            nearest_stations[["msoa11cd", "dist_to_station_m"]], 
            on="msoa11cd", 
            how="left"
        )
        gm_msoas["dist_to_station_km"] = gm_msoas["dist_to_station_m"] / 1000

        # Calculate percentage greenspace
        logger.info("Calculating greenspace percentages...")
        gm_msoas["msoa_area_m2"] = gm_msoas.geometry.area
        
        greenspace_intersections = gpd.overlay(
            gm_msoas[["msoa11cd", "geometry"]], 
            greenspaces[["geometry"]], 
            how="intersection"
        )
        greenspace_intersections["green_area_m2"] = greenspace_intersections.geometry.area
        
        greenspace_totals = greenspace_intersections.groupby("msoa11cd", as_index=False)["green_area_m2"].sum()
        gm_msoas = gm_msoas.merge(greenspace_totals, on="msoa11cd", how="left")
        gm_msoas["green_area_m2"] = gm_msoas["green_area_m2"].fillna(0)
        gm_msoas["pct_greenspace"] = gm_msoas["green_area_m2"] / gm_msoas["msoa_area_m2"]

        # Calculate school count within 1km
        logger.info("Calculating school counts within 1km...")
        msoa_centroids = gm_msoas.set_geometry("centroid")[["msoa11cd", "centroid"]].rename_geometry("geometry")
        
        # Create 1km buffers around centroids
        school_buffers = msoa_centroids.copy()
        school_buffers["geometry"] = school_buffers.geometry.buffer(1000)  # 1km in meters
        
        # Count schools within buffers
        school_intersections = gpd.sjoin(schools[["geometry"]], school_buffers, predicate="within", how="left")
        school_counts = school_intersections.groupby("msoa11cd", as_index=False).size()
        school_counts = school_counts.rename(columns={"size": "school_count_1km"})
        
        gm_msoas = gm_msoas.merge(school_counts, on="msoa11cd", how="left")
        gm_msoas["school_count_1km"] = gm_msoas["school_count_1km"].fillna(0)

        # Export results
        output_columns = ["msoa11cd", "dist_to_station_km", "pct_greenspace", "school_count_1km"]
        gm_msoas[output_columns].to_csv(args.out, index=False)
        logger.info(f"Successfully wrote spatial controls to {args.out}")

    except Exception as e:
        logger.error(f"Error computing spatial controls: {e}")
        raise SpatialDataError(f"Failed to compute spatial controls: {e}")


def compute_msoa_areas(args: argparse.Namespace) -> None:
    """
    Compute MSOA areas in square kilometers.
    
    Args:
        args: Command line arguments with input/output paths
    """
    args.out.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Computing MSOA areas...")

    try:
        gdf = gpd.read_file(args.shapefile).to_crs(epsg=EPSG_BNG)
        valid_codes = set(read_msoa_codes(args.gm_codes))
        
        gm_msoas = gdf[gdf["geo_code"].isin(valid_codes)].copy()
        gm_msoas["MSOA11CD"] = gm_msoas["geo_code"]
        gm_msoas["area_km2"] = gm_msoas.geometry.area / 1_000_000  # Convert m² to km²
        
        result_columns = ["MSOA11CD", "area_km2"]
        gm_msoas[result_columns].to_csv(args.out, index=False)
        logger.info(f"Successfully computed areas for {len(gm_msoas)} MSOAs")

    except Exception as e:
        logger.error(f"Error computing MSOA areas: {e}")
        raise SpatialDataError(f"Failed to compute MSOA areas: {e}")


def compute_cbd_distances(args: argparse.Namespace) -> None:
    """
    Calculate distance from MSOA centroids to Central Business District.
    
    Args:
        args: Command line arguments including CBD coordinates and paths
    """
    try:
        from shapely.geometry import Point
        from pyproj import Transformer
    except ImportError as e:
        raise RuntimeError(f"Required dependencies not available: {e}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Computing distances to CBD at ({args.cbd_lat}, {args.cbd_lon})...")

    try:
        msoas = gpd.read_file(args.shapefile).to_crs(epsg=EPSG_BNG)
        valid_codes = set(read_msoa_codes(args.gm_codes))
        gm_msoas = msoas[msoas["geo_code"].isin(valid_codes)].copy()

        # Transform CBD coordinates to British National Grid
        transformer = Transformer.from_crs(EPSG_WGS84, EPSG_BNG, always_xy=True)
        cbd_x, cbd_y = transformer.transform(args.cbd_lon, args.cbd_lat)
        cbd_point = Point(cbd_x, cbd_y)

        # Calculate distances
        gm_msoas["centroid"] = gm_msoas.geometry.centroid
        gm_msoas["dist_to_CBD_km"] = gm_msoas["centroid"].distance(cbd_point) / 1000.0

        result_columns = ["geo_code", "dist_to_CBD_km"]
        gm_msoas[result_columns].to_csv(args.out, index=False)
        logger.info(f"Successfully computed CBD distances for {len(gm_msoas)} MSOAs")

    except Exception as e:
        logger.error(f"Error computing CBD distances: {e}")
        raise SpatialDataError(f"Failed to compute CBD distances: {e}")


def compute_historical_station_distances(args: argparse.Namespace) -> None:
    """
    Calculate distance to nearest station for multiple years based on opening dates.
    
    Args:
        args: Command line arguments with years, paths, and station opening data
    """
    args.out.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Computing station distances for years: {args.years}")

    try:
        # Load MSOA data
        msoas = gpd.read_file(args.shapefile).to_crs(epsg=EPSG_BNG)
        valid_codes = set(read_msoa_codes(args.gm_codes))
        gm_msoas = msoas[msoas["geo_code"].isin(valid_codes)].copy()
        gm_msoas = gm_msoas.rename(columns={"geo_code": "msoa11cd"})

        # Load station data with opening years
        stations = gpd.read_file(args.stations).to_crs(epsg=EPSG_BNG)
        opening_data = pd.read_csv(args.opening_csv)
        validate_required_columns(opening_data, ["name", "opening_year"])
        
        stations = stations.merge(opening_data, on="name", how="inner")
        stations = stations.dropna(subset=["opening_year"])
        stations["geometry"] = stations.geometry.centroid
        logger.info(f"Loaded {len(stations)} stations with opening year data")

        # Prepare output dataframe
        msoa_geometries = gm_msoas[["msoa11cd", "geometry"]].copy()
        results = gm_msoas[["msoa11cd"]].copy()

        # Calculate distances for each year
        for year in args.years:
            logger.info(f"Processing year {year}...")
            year_stations = stations[stations["opening_year"] <= year]
            logger.info(f"  {len(year_stations)} stations available in {year}")
            
            if len(year_stations) > 0:
                nearest = gpd.sjoin_nearest(
                    msoa_geometries, 
                    year_stations[["geometry"]], 
                    how="left", 
                    distance_col=f"dist_{year}_m"
                )
                results[f"dist_to_station_{year}"] = nearest[f"dist_{year}_m"] / 1000.0
            else:
                logger.warning(f"No stations available for year {year}")
                results[f"dist_to_station_{year}"] = float('inf')

        results.to_csv(args.out, index=False)
        logger.info(f"Successfully computed historical station distances")

    except Exception as e:
        logger.error(f"Error computing historical station distances: {e}")
        raise SpatialDataError(f"Failed to compute historical station distances: {e}")


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure the command line argument parser.
    
    Returns:
        Configured ArgumentParser instance
    """
    base_path = get_base_directory()
    
    parser = argparse.ArgumentParser(
        prog="gm_tools.py",
        description="Greater Manchester spatial data processing toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s extract-osm --north 53.6 --south 53.3
  %(prog)s spatial-controls --out results.csv
  %(prog)s area --shapefile custom_boundaries.shp
        """
    )
    
    subparsers = parser.add_subparsers(dest="command", required=True, help="Available commands")

    # Extract OSM features command
    osm_parser = subparsers.add_parser(
        "extract-osm", 
        help="Extract OpenStreetMap features for a specified bounding box"
    )
    osm_parser.add_argument("--north", type=float, default=53.6, 
                           help="Northern boundary of bounding box (default: 53.6)")
    osm_parser.add_argument("--south", type=float, default=53.3,
                           help="Southern boundary of bounding box (default: 53.3)")
    osm_parser.add_argument("--east", type=float, default=-1.95,
                           help="Eastern boundary of bounding box (default: -1.95)")
    osm_parser.add_argument("--west", type=float, default=-2.60,
                           help="Western boundary of bounding box (default: -2.60)")
    osm_parser.add_argument("--outdir", type=Path, default=base_path/"data/raw/osm",
                           help="Output directory for OSM files")
    osm_parser.set_defaults(func=extract_osm_features)

    # Spatial controls command
    controls_parser = subparsers.add_parser(
        "spatial-controls", 
        help="Compute spatial control variables (distance to station, greenspace %, school count)"
    )
    controls_parser.add_argument("--shapefile", type=Path, 
                                default=base_path/"data/raw/shapefiles/infuse_msoa_lyr_2011.shp",
                                help="Path to MSOA shapefile")
    controls_parser.add_argument("--gm-codes", type=Path, 
                                default=base_path/"data/reference/msoa_codes.txt",
                                help="Path to file containing MSOA codes")
    controls_parser.add_argument("--stations", type=Path, 
                                default=base_path/"data/raw/osm/transport_stations.geojson",
                                help="Path to stations GeoJSON file")
    controls_parser.add_argument("--greenspace", type=Path, 
                                default=base_path/"data/raw/osm/greenspaces.geojson",
                                help="Path to greenspaces GeoJSON file")
    controls_parser.add_argument("--schools", type=Path, 
                                default=base_path/"data/raw/osm/schools.geojson",
                                help="Path to schools GeoJSON file")
    controls_parser.add_argument("--out", type=Path, 
                                default=base_path/"data/processed/spatial/spatial_controls.csv",
                                help="Output CSV file path")
    controls_parser.set_defaults(func=compute_spatial_controls)

    # Area computation command
    area_parser = subparsers.add_parser("area", help="Compute MSOA areas in square kilometers")
    area_parser.add_argument("--shapefile", type=Path, 
                            default=base_path/"data/raw/shapefiles/infuse_msoa_lyr_2011.shp",
                            help="Path to MSOA shapefile")
    area_parser.add_argument("--gm-codes", type=Path, 
                            default=base_path/"data/reference/msoa_codes.txt",
                            help="Path to file containing MSOA codes")
    area_parser.add_argument("--out", type=Path, 
                            default=base_path/"data/processed/spatial/msoa_areas_km2.csv",
                            help="Output CSV file path")
    area_parser.set_defaults(func=compute_msoa_areas)

    # CBD distance command
    cbd_parser = subparsers.add_parser("dist-cbd", help="Calculate distance from MSOA centroids to CBD")
    cbd_parser.add_argument("--shapefile", type=Path, 
                           default=base_path/"data/raw/shapefiles/infuse_msoa_lyr_2011.shp",
                           help="Path to MSOA shapefile")
    cbd_parser.add_argument("--gm-codes", type=Path, 
                           default=base_path/"data/reference/msoa_codes.txt",
                           help="Path to file containing MSOA codes")
    cbd_parser.add_argument("--cbd-lat", type=float, default=53.4794,
                           help="CBD latitude coordinate (default: 53.4794)")
    cbd_parser.add_argument("--cbd-lon", type=float, default=-2.2453,
                           help="CBD longitude coordinate (default: -2.2453)")
    cbd_parser.add_argument("--out", type=Path, 
                           default=base_path/"data/processed/spatial/msoa_cbd_distances.csv",
                           help="Output CSV file path")
    cbd_parser.set_defaults(func=compute_cbd_distances)

    # Historical station distances command
    historical_parser = subparsers.add_parser(
        "dist-station-over-time", 
        help="Calculate nearest station distances for multiple years"
    )
    historical_parser.add_argument("--years", nargs="+", type=int, default=[1991, 2001, 2011, 2021],
                                  help="Years to analyze (default: 1991 2001 2011 2021)")
    historical_parser.add_argument("--shapefile", type=Path, 
                                  default=base_path/"data/raw/shapefiles/infuse_msoa_lyr_2011.shp",
                                  help="Path to MSOA shapefile")
    historical_parser.add_argument("--gm-codes", type=Path, 
                                  default=base_path/"data/reference/msoa_codes.txt",
                                  help="Path to file containing MSOA codes")
    historical_parser.add_argument("--stations", type=Path, 
                                  default=base_path/"data/raw/osm/transport_stations.geojson",
                                  help="Path to stations GeoJSON file")
    historical_parser.add_argument("--opening-csv", type=Path, 
                                  default=base_path/"data/reference/station_opening_years.csv",
                                  help="Path to CSV with station opening years")
    historical_parser.add_argument("--out", type=Path, 
                                  default=base_path/"data/processed/spatial/historical_station_distances.csv",
                                  help="Output CSV file path")
    historical_parser.set_defaults(func=compute_historical_station_distances)

    return parser


def main() -> None:
    """Main entry point for the CLI application."""
    try:
        parser = create_argument_parser()
        args = parser.parse_args()
        
        logger.info(f"Starting {args.command} command...")
        args.func(args)
        logger.info(f"Command {args.command} completed successfully")
        
    except KeyboardInterrupt:
        logger.info("Operation cancelled by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
