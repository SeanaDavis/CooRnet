#' @import igraph
#' @importFrom dplyr group_by mutate %>% summarize
#' @import longurl

build_coord_graph <- function(ct_shares.df, coordinated_shares, percentile_edge_weight=0.90, timestamps=FALSE) {

  ###########################################
  # Get coordinated entities and network ####
  ###########################################

  cat("\nBuilding the graph...")

  el <- coordinated_shares[,c(3,5,4)] # drop unnecesary columns
  el$account.url <- trimws(el$account.url) # remove white space from platform.id
  v1 <- data.frame(node=unique(el$account.url), type=1) # create a dataframe with nodes and type 0=url 1=page
  v2 <- data.frame(node=unique(el$url), type=0)
  v <- rbind(v1,v2)

  g2.bp <- graph.data.frame(el,directed = T, vertices = v) # makes the biap
  g2.bp <- igraph::simplify(g2.bp, remove.multiple = T, remove.loops = T ,edge.attr.comb = "min") # simply the bipartite netwrok to avoid problems with resulting edge weight in projected network
  full_g <- suppressWarnings(bipartite.projection(g2.bp, multiplicity = T)$proj2) # project page-page network

  all_account_info <- ct_shares.df %>%
    group_by(account.url) %>%
    summarize(shares = n(),
              avg.account.subscriberCount=mean(account.subscriberCount))

  # group the pages that changed names or handles
  ct_shares.df <- ct_shares.df %>%
    group_by(account.url) %>%
    mutate(name.changed = ifelse(length(unique(account.name))>1, TRUE, FALSE),
           handle.changed = ifelse(length(unique(account.handle))>1, TRUE, FALSE),
           account.name = paste(unique(account.name), collapse = " | "),
           account.handle = paste(unique(account.handle), collapse = " | "))

  more.account.info <- ct_shares.df[, c("account.id", "account.name", "name.changed", "handle.changed", "account.handle",
                                        "account.url", "account.platform", "account.platformId", "account.verified")]

  rm(ct_shares.df, coordinated_shares)

  more.account.info <- unique(more.account.info)
  all_account_info <- merge(all_account_info, more.account.info, by="account.url")

  # add vertex attributes
  vertex.info <- subset(all_account_info, as.character(all_account_info$account.url) %in% V(full_g)$name)

  V(full_g)$shares <- sapply(V(full_g)$name, function(x) vertex.info$shares[vertex.info$account.url == x])
  V(full_g)$avg.account.subscriberCount <- sapply(V(full_g)$name, function(x) vertex.info$avg.account.subscriberCount[vertex.info$account.url == x])
  V(full_g)$account.platform <- sapply(V(full_g)$name, function(x) vertex.info$account.platform[vertex.info$account.url == x])
  V(full_g)$account.name <- sapply(V(full_g)$name, function(x) vertex.info$account.name[vertex.info$account.url == x])
  V(full_g)$account.verified <- sapply(V(full_g)$name, function(x) vertex.info$account.verified[vertex.info$account.url == x])
  V(full_g)$account.handle <- sapply(V(full_g)$name, function(x) vertex.info$account.handle[vertex.info$account.url == x])
  V(full_g)$name.changed <- sapply(V(full_g)$name, function(x) vertex.info$name.changed[vertex.info$account.url == x])
  V(full_g)$handle.changed <- sapply(V(full_g)$name, function(x) vertex.info$handle.changed[vertex.info$account.url == x])

  # keep only highly coordinated entities
  V(full_g)$degree <- degree(full_g)
  q <- quantile(E(full_g)$weight, percentile_edge_weight) # set the percentile_edge_weight number of repetedly coordinated link sharing to keep
  highly_connected_g <- induced_subgraph(graph = full_g, vids = V(full_g)[V(full_g)$degree > 0 ]) # filter for degree
  highly_connected_g <- subgraph.edges(highly_connected_g, eids = which(E(highly_connected_g)$weight >= q),delete.vertices = T) # filter for edge weight

  if (timestamps==TRUE) {
    cat("\n\nAdding timestamps. Please be patient... :)")

    EL <- as.data.frame(as_edgelist(highly_connected_g))
    EL$weight <- E(highly_connected_g)$weight
    EL$coord_shares <- 0
    names(EL) <- c("V1","V2","weight","coord_share")

    coord_shares <- list()

    EL <- as_tibble(EL)

    ts_on_edge <- function(x,output){
      shared_2 <- intersect(neighbors(graph = g2.bp,v = V(g2.bp)[name==x[1]],mode = "out"),
                            neighbors(graph = g2.bp,v = V(g2.bp)[name==x[2]],mode = "out"))
      cs <- rep(NA,length(shared_2))
      cs <- sapply(shared_2, function(i) E(g2.bp)[get.edge.ids(graph = g2.bp,directed = F,vp = c(i,V(g2.bp)[name==x[1]]))]$share_date)

      return(cs)
    }

    time1 <- Sys.time()
    E(highly_connected_g)$t_coord_share <- apply(EL,1,ts_on_edge)
    time2 <- Sys.time()
    time2-time1
  }

  # find and annotate nodes-components
  V(highly_connected_g)$component <- components(highly_connected_g)$membership
  V(highly_connected_g)$degree <- degree(highly_connected_g) # re-calculate the degree on the subgraph
  V(highly_connected_g)$strength <- strength(highly_connected_g) # sum up the edge weights of the adjacent edges for each vertex

  highly_connected_coordinated_entities <- igraph::as_data_frame(highly_connected_g, "vertices")

  # detects entities type
  cat("\nDetecting entities type...")

  highly_connected_coordinated_entities <- highly_connected_coordinated_entities %>%
    mutate(account.url_expanded = ifelse(grepl("/groups/", name), name, longurl::expand_urls(name, warn = F, agent = "cooRnet-r-package", .progress = TRUE)$expanded_url),
           entity_type = ifelse(grepl("/groups/", account.url_expanded), "group", "page")) %>%
    select(-account.url_expanded)

  # copy type into graph
  V(highly_connected_g)$entity_type <- sapply(V(highly_connected_g)$name, function(x) highly_connected_coordinated_entities$entity_type[highly_connected_coordinated_entities$name == x])

  rownames(highly_connected_coordinated_entities) <- 1:nrow(highly_connected_coordinated_entities)
  colnames(more.account.info)[5] <- "name" # please use column name and not number

  highly_c_list <- list(highly_connected_g, highly_connected_coordinated_entities, q)

  cat("\nDone!")

  return(highly_c_list)
}
