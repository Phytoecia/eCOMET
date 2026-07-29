# ExportCytoscape

Export a `mmo` as a node table and edge table for visualisation in
Cytoscape. Edges are derived from a stored feature dissimilarity matrix.
Nodes carry abundance statistics and any annotation columns present in
`mmo$sirius_annot` and `mmo$feature_info`.

## Usage

``` r
ExportCytoscape(
  mmo,
  distance = "dreams",
  outprefix = "cytoscape_export",
  sim_threshold = 0.7,
  top_k = NULL,
  group_col = "group",
  sample_col = "sample"
)
```

## Arguments

  - mmo:
    
    `mmo`. Must contain `mmo$feature_data` and `mmo$metadata`.
    `mmo$feature_info` and `mmo$sirius_annot` are used when present.

  - distance:
    
    Name of the dissimilarity matrix to use for edges (default:
    `"dreams"`). Passed to `GetDistanceMat()`.

  - outprefix:
    
    File path prefix for output CSVs (default: `"cytoscape_export"`).

  - sim\_threshold:
    
    Minimum pairwise similarity (= 1 - dissimilarity) for an edge to be
    retained (default: `0.7`). Set to `0` to disable.

  - top\_k:
    
    Integer or `NULL`. If not `NULL`, each node retains only its `top_k`
    most similar neighbours after threshold filtering (default: `NULL` =
    no k limit).

  - group\_col:
    
    Metadata column defining groups for per-group mean abundance columns
    (default: `"group"`).

  - sample\_col:
    
    Metadata column matching sample IDs to feature columns (default:
    `"sample"`).

## Value

Invisibly returns a named list: `nodes` (data.frame), `edges`
(data.frame), `nodes_path`, `edges_path` (file paths).

## Details

**Edge filtering:** Two complementary filters control which edges are
retained. Both are applied together (a pair must pass both to be
included):

  - `sim_threshold` – drop all edges with similarity below this value.
    Higher values produce sparser, higher-confidence networks.

  - `top_k` – for each node, retain only its `top_k` most similar
    neighbours (by similarity). `NULL` keeps all edges that pass the
    threshold. This is useful for preventing highly-connected hub nodes
    from dominating the layout.

A warning is issued when the retained edge count exceeds 50 000, as
large networks can be slow to render in Cytoscape.

**Node table columns:** The node table always includes `id`,
`prevalence` (proportion of samples with detected abundance), and one
`mean_<group>` column per group. Any columns present in
`mmo$feature_info` and `mmo$sirius_annot` are appended automatically –
the function does not assume specific column names. Column names are
sanitised to replace characters that cause problems in Cytoscape
(spaces, `#`, `:`).

**Loading in Cytoscape:**

1.  File -\> Import -\> Network from File -\> select `_edges.csv`. Map
    `source` as Source Node, `target` as Target Node, `similarity` as
    Edge Attribute.

2.  File -\> Import -\> Table from File -\> select `_nodes.csv`. Map
    `id` as the Key column matching node names.

3.  In the Style panel, map `Fill Color` to any annotation column (e.g.
    `NPC_pathway`) to color nodes by compound class.

## Examples

``` r
if (FALSE) {
# Default: similarity >= 0.7
ExportCytoscape(mmo, distance = "dreams",
                outprefix = "output/network")

# Sparser network: top 5 neighbours per node, similarity >= 0.6
ExportCytoscape(mmo, distance = "dreams",
                outprefix = "output/network_k5",
                sim_threshold = 0.6, top_k = 5)
}
```
